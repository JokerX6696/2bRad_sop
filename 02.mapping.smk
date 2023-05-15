rule mkref_all_fa:
    input:
        data_list = 'result/01.QC/Prepare/ustacks.list'
    output:
        ref = 'result/02.Mapping/mkref/all.fa.gz',
    params:
        
    log:
        o = "logs/02.Mapping/mkref_all_fa/o.mkref_all_fa.log",
        e = "logs/02.Mapping/mkref_all_fa/e.mkref_all_fa.log"
    benchmark: 
        "benchmarks/02.Mapping/mkref_all_fa/mkref_all_fa.txt"
    #threads: cores
    shell:
        """
        {perl} \
        {mkref_prepare} \
        -l {input.data_list} \
        -o result/02.Mapping/mkref \
        1>{log.o} 2>{log.e}
        """

rule mkref_ustacks:
    input:
        fa = rules.mkref_all_fa.output.ref
    output:
        ref_tag = 'result/02.Mapping/mkref/all.tags.tsv.gz'       
    params:
        
    log:
        o = "logs/02.Mapping/mkref/o.mkref_ustacks.log",
        e = "logs/02.Mapping/mkref/e.mkref_ustacks.log"
    benchmark: 
        "benchmarks/02.Mapping/mkref/mkref_ustacks.txt"
    #threads: cores
    shell:
        """
        {ustack} \
        -m 3 -M 2 -t gzfasta -H -p 40 \
        -f {input.fa} \
        -o result/02.Mapping/mkref \
        1>{log.o} 2>{log.e}
        """
        
rule mkref_ref_HQ_fa:
    input:
        tags = rules.mkref_ustacks.output.ref_tag
    output: 
        HQ_fa = 'result/02.Mapping/mkref/ref_HQ.fa.gz'
  
    params:
        
    log:
        o = "logs/02.Mapping/mkref/o.mkref_ref_HQ_fa.log",
        e = "logs/02.Mapping/mkref/e.mkref_ref_HQ_fa.log"
    benchmark: 
        "benchmarks/02.Mapping/mkref/mkref.txt"
    #threads: cores
    shell:
        """
        {perl} \
        {mkref_stat_and_deal} \
        -i {input.tags} \
        -o result/02.Mapping/mkref \
        1>>{log.o} 2>>{log.e}
        """

####################################### 有参 or 无参 ref 调整 1 ##################################################
if refed == 'no':
    fa = rules.mkref_ref_HQ_fa.output.HQ_fa
else:
    fa = ref_fa

#################################################################################################################

rule map_enzyme_index:
    input:
        ref_HQ = fa
    output:
        fa = 'result/02.Mapping/map_enzyme_index/ref.fa',
        idx = "result/02.Mapping/map_enzyme_index/ref.fa.index.sai"
    params:
        prefix = 'result/02.Mapping/map_enzyme_index/ref'
        
    log:
        o = "logs/02.Mapping/map_enzyme_index/o.map_enzyme_index.log",
        e = "logs/02.Mapping/map_enzyme_index/e.map_enzyme_index.log"
    benchmark: 
        "benchmarks/02.Mapping/map_enzyme_index/map_enzyme_index.txt"
    #threads: cores
    shell:
        """
        {perl} {EnzymeSite} -g {input.ref_HQ} -s 3 -p {params.prefix} 1>{log.o} 2>{log.e}
        {builder} {output.fa} 1>>{log.o} 2>>{log.e}
        """

rule unzip_fq:
    input:
        fq = "clean_data/{sample}." + enzyme + ".fq.gz"
    output:
        "result/02.Mapping/unzip/{sample}.fastq"
    params:

    log:
        o = "logs/02.Mapping/unzip/o.{sample}_map_soap.log",
        e = "logs/02.Mapping/unzip/e.{sample}_map_soap.log"
    benchmark: 
        "benchmarks/02.Mapping/unzip/{sample}_unzip.txt"
    #threads: cores
    shell:  # 这里 1 本身就是 解压后的输出 所以不能将 1 重定向到 log.o
        """
        gzip -dc {input.fq} > {output} 2>{log.e}
        
        """


rule map_soap:
    input:
        fq = rules.unzip_fq.output
    output:
        map_tmp = temp("result/02.Mapping/soap/{sample}.map_tmp.gz")
    params:
        
    log:
        o = "logs/02.Mapping/soap/o.{sample}_map_soap.log",
        e = "logs/02.Mapping/soap/{sample}.log"
    benchmark: 
        "benchmarks/02.Mapping/soap/{sample}_map_soap.txt"
    #threads: cores
    shell:
        """
        /data/software/install/soap2.2.21/soap \
        -a {input.fq} \
        -D "result/02.Mapping/map_enzyme_index/ref.fa.index" \
        -M 4 -r 0 -v 2 \
        -o {output} \
        -p 10 1>{log.o} 2>{log.e}
        """

rule map_filter:
    input:
        maped = rules.map_soap.output.map_tmp
    output:
        temp("result/02.Mapping/soap/{sample}.soap")
    params:
        
    log:
        o = "logs/02.Mapping/soap/o.{sample}_map_filter.log",
        e = "logs/02.Mapping/soap/e.{sample}_map_filter.log"
    benchmark: 
        "benchmarks/02.Mapping/soap/{sample}_map_filter.txt"
    #threads: cores
    shell:
        """
        awk '$9==1' {input.maped} > {output} 2>{log.e}
        """

rule readsmap:
    input:
        rules.map_filter.output
    output:
        "result/02.Mapping/soap/{sample}.map.gz"
    params:
        
    log:
        o = "logs/02.Mapping/soap/o.{sample}_readsmap.log",
        e = "logs/02.Mapping/soap/e.{sample}_readsmap.log"
    benchmark: 
        "benchmarks/02.Mapping/soap/{sample}_readsmap.txt"
    #threads: cores
    shell:
        """
        /data/software/install/miniconda3/envs/perl.5.26.2/bin/perl \
        /data/scripts/rad/2brad/01.pipline/version4.1/script/readsmap.pl \
        -i {input} \
        -r {fa} \
        -o {output} 2>{log.e}

        """

rule map_stat:
    input:
        expand("logs/02.Mapping/soap/{sample}.log",sample=samples)
    output:
        "result/02.Mapping/map_stat/soap_rate.xls"
    params:
        
    log:
        o = "logs/02.Mapping/map_stat/o.map_stat.log",
        e = "logs/02.Mapping/map_stat/e.map_stat.log"
    benchmark: 
        "benchmarks/02.Mapping/map_stat/map_stat.txt"
    #threads: cores
    shell:
        """
        /data/software/install/miniconda3/envs/perl.5.26.2/bin/perl \
        /data/scripts/rad/2brad/01.pipline/version4.1/script/soap_ratio_stat.pl \
        -id {input} \
        -o {output}
        """
   
