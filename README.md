# **血浆微生物reads中鉴定结构保守的ncRNA**

![image-20251229202952292](/Users/qin/Library/Application Support/typora-user-images/image-20251229202952292.png)

### 数据前处理（去除接头和低质量reads）+ Mapping

跑dt2-BC6-pe.snakefile

```
######################################################################

#config
#configfile: "/lulabdata/zhanqing/projects/SLE_primary_disease/output-PE/script-PE/cnode-config.yaml"

#reference
tool_dir="/lulabdata/zhanqing/BioII/cfRNAseq/tools"
index_dir="/lulabdata/zhanqing/BioII/cfRNAseq/reference/star"
fasta_dir="/lulabdata/zhanqing/BioII/cfRNAseq/reference/fasta"
bed_dir="/lulabdata/zhanqing/BioII/cfRNAseq/reference/bed/gencode.v38.sorted"
gtf_dir="/lulabdata/zhanqing/BioII/cfRNAseq/reference/gtf"
gene_table_dir="/lulabdata/zhanqing/BioII/cfRNAseq/reference/gene-table"
bt2_dir="/lulabdata/zhanqing/BioII/cfRNAseq/reference/bowtie2_index"


input_dir="/lulabdata2/SLE_primary_disease/output"
output_dir="/lulabdata3/liuke/SLE/output-PE"
sample_id_file="/lulabdata/liuke/pipeline/SLE/SLE_Primary_Disease/mapping/sample_rest.txt"

#sample_id_file
file=open(expand('{sample_id_file}',sample_id_file=sample_id_file)[0])
sample_ids=[]
for line in file:
	sample_ids.append(line.strip())
sample_id_1=sample_ids[0]

#fastqc
thread_fastqc=12
#cutadapt
thread_cutadapt=16
clean_levels=['cutadapt','trimGC']
R1_3adapter="NNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
R1_5adapter="T{150}VN"
R2_3adapter="TGCTNNNNNNNNATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
R2_3adapter_polyA="AAAAAAAAAA"

#trimGC
strandness="reverse"


#mapping
thread_mapping=16
threads_decompress=16
#map_steps=['spikein_long','univec','rRNA','hg38_v38','circRNA']
map_steps=['spikein_long','univec','rRNA','hg38_v38']
map_steps_sortbyCoord=['hg38_v38']
map_steps_dedup=['hg38_v38']
map_steps_sortbyName=['hg38_v38']
#kraken
taxoLevel: ['P','G','D']
#count
count_levels=['gencode']
#genome_region
priority=['MT_rRNA,MT_mRNA,MT_tRNA,chrM.all,mRNA,lncRNA,pseudogene,tRNA,srpRNA,snoRNA,snRNA,Y_RNA,misc_RNA,tucpRNA,exon,intron,intergenic']
######################################################################
######################################################################

#pipeline
def get_all_inputs(wildcards):
	available_inputs = dict(
		cutadapt=expand('{input_dir}/cutadapt/{sample_id}_{mate_index}.fastq.gz',
                        input_dir=input_dir, sample_id=sample_ids, mate_index=[1, 2]),
        trimGC=expand('{output_dir}/trimGC/{sample_id}_{mate_index}.fastq.gz',
                        output_dir=output_dir, sample_id=sample_ids, mate_index=[1, 2]),
		fastqc_clean=expand('{output_dir}/fastqc_clean/{sample_id}_{mate_index}_fastqc.{type}',
                        output_dir=output_dir, sample_id=sample_ids, mate_index=[1, 2], type=['html','zip']),
		bam=expand('{output_dir}/bam/{sample_id}/{map_step}.bam',
                        output_dir=output_dir, sample_id=sample_ids, map_step=map_steps)
	)
	enabled_inputs = list(available_inputs.keys())
	inputs = []
	for key, l in available_inputs.items():
		if key in enabled_inputs:
			inputs += l
	return inputs


rule all:
	input:
		get_all_inputs



rule trimGC_pe:
		input:
				fastq1=input_dir+'/cutadapt/{sample_id}_1.fastq.gz',
				fastq2=input_dir+'/cutadapt/{sample_id}_2.fastq.gz'
		output:
				fastq1='{output_dir}/trimGC/{sample_id}_1.fastq.gz',
				fastq2='{output_dir}/trimGC/{sample_id}_2.fastq.gz'
		log:
				'{output_dir}/log/{sample_id}/trimGC.log'
		shell:
				'''
				python {tool_dir}/trimGC_zq.py -m 16 -s {strandness} \
                -o {output_dir}/trimGC/{wildcards.sample_id} \
                -i {input_dir}/cutadapt/{wildcards.sample_id} > {log} 2>&1
				'''

rule fastqc_clean:
	input:
		fastq1='{output_dir}/trimGC/{sample_id}_1.fastq.gz',
		fastq2='{output_dir}/trimGC/{sample_id}_2.fastq.gz'
	output:
		html1='{output_dir}/fastqc_clean/{sample_id}_1_fastqc.html',
		html2='{output_dir}/fastqc_clean/{sample_id}_2_fastqc.html',
		zip1='{output_dir}/fastqc_clean/{sample_id}_1_fastqc.zip',
		zip2='{output_dir}/fastqc_clean/{sample_id}_2_fastqc.zip'
	log:
		'{output_dir}/log/{sample_id}/fastqc_clean.log'
	shell:
		'''
		fastqc --threads {thread_fastqc} -o {output_dir}/fastqc_clean {input.fastq1} > {log} 2>&1
		fastqc --threads {thread_fastqc} -o {output_dir}/fastqc_clean {input.fastq2} > {log} 2>&1
		'''

map_command_pe = '''STAR --genomeDir {params.index} \
              --readFilesIn {input.fastq1} {input.fastq2} \
              --runThreadN {thread_mapping} \
              --outFileNamePrefix {params.output_prefix} \
              --outSAMtype BAM Unsorted \
              --outReadsUnmapped Fastx \
              --readFilesCommand gzip -d -c \
              --outSAMmultNmax -1 \
	      --seedPerWindowNmax {params.seedPerWindowNmax} > {log} 2>&1
	  mv {params.output_prefix}Aligned.out.bam {output.bam}
	  {tool_dir}/bbmap/repair.sh -Xms1000m -Xmx100000m overwrite=t \
		in={params.output_prefix}Unmapped.out.mate1 \
		in2={params.output_prefix}Unmapped.out.mate2 \
		out={params.output_prefix}Unmapped.out.mate3 \
		out2={params.output_prefix}Unmapped.out.mate4 >> {log} 2>&1
	  pigz -c {params.output_prefix}Unmapped.out.mate3 > {output.unmapped1}
	  pigz -c {params.output_prefix}Unmapped.out.mate4 > {output.unmapped2}
	  rm -f {params.output_prefix}Unmapped.out.mate1 {params.output_prefix}Unmapped.out.mate2 {params.output_prefix}Unmapped.out.mate3 {params.output_prefix}Unmapped.out.mate4
        '''

map_command_pe_mismatch = '''STAR --genomeDir {params.index} \
              --readFilesIn {input.fastq1} {input.fastq2} \
              --runThreadN {thread_mapping} \
              --outFileNamePrefix {params.output_prefix} \
              --outSAMtype BAM Unsorted \
              --outReadsUnmapped Fastx \
              --twopassMode Basic \
              --readFilesCommand gzip -d -c \
              --outSAMmultNmax -1 \
              --outSAMstrandField intronMotif \
              --outSAMattributes All \
              --outFilterType BySJout\
              --outFilterMismatchNoverLmax 0.04 \
              --alignIntronMax 1000000 \
              --alignMatesGapMax 1000000 \
	      --seedPerWindowNmax {params.seedPerWindowNmax} > {log} 2>&1
	  mv {params.output_prefix}Aligned.out.bam {output.bam}
	  {tool_dir}/bbmap/repair.sh -Xms1000m -Xmx100000m overwrite=t \
		in={params.output_prefix}Unmapped.out.mate1 \
		in2={params.output_prefix}Unmapped.out.mate2 \
		out={params.output_prefix}Unmapped.out.mate3 \
		out2={params.output_prefix}Unmapped.out.mate4 >> {log} 2>&1
	  pigz -c {params.output_prefix}Unmapped.out.mate3 > {output.unmapped1}
	  pigz -c {params.output_prefix}Unmapped.out.mate4 > {output.unmapped2}
	  rm -f {params.output_prefix}Unmapped.out.mate1 {params.output_prefix}Unmapped.out.mate2 {params.output_prefix}Unmapped.out.mate3 {params.output_prefix}Unmapped.out.mate4
        '''

#-Xms1000m -Xmx100000m !!!!!!!!
#!!!!!!   --outSAMmultNmax -1   !!!!!!


rule map_spikein_long:
	input:
		fastq1='{output_dir}/trimGC/{sample_id}_1.fastq.gz',
		fastq2='{output_dir}/trimGC/{sample_id}_2.fastq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/spikein_long.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/spikein_long_1.fq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/spikein_long_2.fq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star-pe/spikein_long/mapping.log'
	params:
		index=index_dir+'/spikein_long',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star-pe/spikein_long/',
		seedPerWindowNmax=20
	run:
		shell(map_command_pe)

rule map_univec:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/spikein_long_1.fq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/spikein_long_2.fq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/univec.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/univec_1.fq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/univec_2.fq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/univec/mapping.log'
	params:
		index=index_dir+'/univec',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/univec/',
		seedPerWindowNmax=20
	run:
	    shell(map_command_pe)

rule map_rRNA:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/univec_1.fq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/univec_2.fq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/rRNA.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/rRNA_1.fq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/rRNA_2.fq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/rRNA/mapping.log'
	params:
		index=index_dir+'/rRNA',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/rRNA/',
		seedPerWindowNmax=20
	run:
		shell(map_command_pe)


rule map_hg38_v38:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/rRNA_1.fq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/rRNA_2.fq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/hg38_v38.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/hg38_v38_1.fq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/hg38_v38_2.fq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/hg38_v38/mapping.log'
	params:
		index=index_dir+'/hg38_v38',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/hg38_v38/',
		seedPerWindowNmax=50
	run:
		shell(map_command_pe_mismatch)


```

### Unmapped reads组装contigs

```
sbatch spades_assemble.slurm
```

spades_assemble.slurm

```
#!/bin/bash
#SBATCH --job-name=spades_assemble        # job 名字
#SBATCH --output=spades_%j.out            # 标准输出 (%j = job ID)
#SBATCH --error=spades_%j.err             # 错误输出
#SBATCH --time=96:00:00                   # 最大运行时间，根据实际修改
#SBATCH --cpus-per-task=12                 # 使用的 CPU 核心数（与脚本里面 -t 参数一致）
#SBATCH -p Lulab


source ~/.bashrc
conda activate bio


echo "Running sample assembly ..."
python /lulabdata3/zhangqin/SLE_microbial_ncRNA/spades_assembled_contigs/spades_assemble.py  /lulabdata3/zhangqin/SLE_microbial_ncRNA/spades_assembled_contigs/ALL_sample.txt

echo "Done."
```

spades_assemble.py

```
import os
import sys
from datetime import datetime

sample_id = sys.argv[1]

target_sample = open(sample_id).read().strip().split("\n")
print(target_sample)

file_path = "/lulabdata3/liuke/SLE/output-PE/unmapped"
os.system('mkdir -p assembled_contigs')

for eles in target_sample :
    input_1 = file_path + "/" + eles.strip() + "/hg38_v38_1.fq.gz"
    input_2 = file_path + "/" + eles.strip() + "/hg38_v38_2.fq.gz"
    output = eles.strip() + "_metatrans_assemble_result"
    print(input_1 + "---" + input_2 + "---" + output)
    # 获取当前时间
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    print(current_time )
    # 写入运行日志
    with open('run.log', 'a') as log_file:
        log_file.write(f"Starting rnaspades.py for {eles.strip()} at {current_time}\n\n")
    os.system('rnaspades.py -1 %s -2 %s -o %s -t 12  --memory 200'%(input_1,input_2,output))
    output_files = eles.strip() + "_transcripts.fasta"
    os.system('sed \'s/NODE_/%s_/g\'  %s/transcripts.fasta | sed \'s/_l/-l/g\' | sed \'s/_c/-c/g\' | sed \'s/_g/-g/g\' > %s '%(eles,output,output_files))
    os.system('rm -rf %s  && mv %s  assembled_contigs'%(output,output_files))
```

```
cat assembled_contigs/*fasta  >  combined_transcripts.fasta

samtools faidx combined_transcripts.fasta

cat combined_transcripts.fasta.fai |  cut -f 2 | awk '{print NR"\t"$0}'  | sed   '1i Name\tLength'  >  combined_transcripts.fasta_length.txt
```

### prodigal_prediction

```
sbatch prodigal_prediction.slurm
```

prodigal_prediction.slurm

```
#!/bin/bash
#SBATCH --job-name=prodigal_prediction        # job 名字
#SBATCH --output=%j.out            # 标准输出 (%j = job ID)
#SBATCH --error=%j.err             # 错误输出
#SBATCH --time=96:00:00                   # 最大运行时间，根据实际修改
#SBATCH --cpus-per-task=8                 # 使用的 CPU 核心数（与脚本里面 -t 参数一致）
#SBATCH -p Lulab


source ~/.bashrc
conda activate bio


prodigal -p meta -i  /lulabdata3/zhangqin/SLE_microbial_ncRNA/spades_assembled_contigs/combined_3060_SLE_transcripts.fasta   -a  SLE_prodigal_predicted_result.aa   -d  SLE_prodigal_predicted_result.fasta    -o  SLE_prodigal_predicted_result.gff  -q  -f gff  -c

```

```
grep ">" SLE_prodigal_predicted_result.fasta | cut -d " " -f 1 | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}' | sort | uniq | sort  >  prodigal_predicted_proteins_contigid.txt

grep ">" combined_3060_SLE_transcripts.fasta  | sort | uniq | sort > all_contigid.txt

comm -3 all_contigid.txt  prodigal_predicted_proteins_contigid.txt  > prodigal_unpredicted_proteins_contigid.txt

sbatch seqkit.slurm
```

```
#!/bin/bash
#SBATCH --job-name=seqkit        # job 名字
#SBATCH --output=%j.out            # 标准输出 (%j = job ID)
#SBATCH --error=%j.err             # 错误输出
#SBATCH --time=96:00:00                   # 最大运行时间，根据实际修改
#SBATCH --cpus-per-task=8                 # 使用的 CPU 核心数（与脚本里面 -t 参数一致）
#SBATCH -p Lulab


source ~/.bashrc
conda activate bio

cat prodigal_unpredicted_proteins_contigid.txt  | cut -d ">" -f 2 | seqkit grep -f - ../spades_assembled_contigs/combined_3060_SLE_transcripts.fasta   -j 8  >  prodigal_unannotated_contigs.fasta

```

### **MetaGeneMark_prediction**

```
sbatch MetaGeneMark.slurm
```

```
#!/bin/bash
#SBATCH --job-name=MegaGeneMark        # job 名字
#SBATCH --output=%j.out            # 标准输出 (%j = job ID)
#SBATCH --error=%j.err             # 错误输出
#SBATCH --time=96:00:00                   # 最大运行时间，根据实际修改
#SBATCH --cpus-per-task=12                 # 使用的 CPU 核心数（与脚本里面 -t 参数一致）
#SBATCH -p Lulab


source ~/.bashrc
conda activate bio

gmhmmp -m ~/software/MetaGeneMark/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod -A MetaGeneMark_predicted_genes.fasta -D MetaGeneMark_predicted_proteins.fasta -f G -o MetaGeneMark_predicted_genes.gff /lulabdata3/zhangqin/SLE_microbial_ncRNA/spades_assembled_contigs/combined_3060_SLE_transcripts.fasta
```

```
grep ">" combined_transcripts_with_diseaseInfo.fasta  | sort | uniq | sort > all_contigid.txt

grep ">" MetaGeneMark_predicted_genes.fasta | cut -d ">" -f 3 | sort | uniq > MetaGeneMark_predicted_proteins_contigid.txt

comm -3 all_contigid.txt  MetaGeneMark_predicted_proteins_contigid.txt  > MetaGeneMark_unpredicted_proteins_contigid.txt

sbatch seqkit.slurm
```

```
#!/bin/bash
#SBATCH --job-name=seqkit        # job 名字
#SBATCH --output=%j.out            # 标准输出 (%j = job ID)
#SBATCH --error=%j.err             # 错误输出
#SBATCH --time=96:00:00                   # 最大运行时间，根据实际修改
#SBATCH --cpus-per-task=8                 # 使用的 CPU 核心数（与脚本里面 -t 参数一致）
#SBATCH -p Lulab


source ~/.bashrc
conda activate bio

cat MetaGeneMark_unpredicted_proteins_contigid.txt  | cut -d ">" -f 2 | seqkit grep -f - ../spades_assembled_contigs/combined_3060_SLE_transcripts.fasta   -j 8  >  MetaGeneMark_unpredicted_proteins_contigid.fasta
```

### **Prodigal_MetageneMark_combined_analysis**

```
cp ../prodigal_prediction/prodigal_predicted_proteins_contigid.txt  ../prodigal_prediction/prodigal_unpredicted_proteins_contigid.txt  ./

cp ../MetaGeneMark_prediction/MetaGeneMark_predicted_proteins_contigid.txt  ../MetaGeneMark_prediction/MetaGeneMark_unpredicted_proteins_contigid.txt  ./

cp  ../MetaGeneMark_prediction/all_contigid.txt  ./

cp combined_transcripts.fasta ./

cat MetaGeneMark_predicted_proteins_contigid.txt  prodigal_predicted_proteins_contigid.txt  | sort | uniq  >  01_MetaGeneMark_Prodigal_predicted_proteins_contigid.txt

seqkit grep -f 01_MetaGeneMark_Prodigal_predicted_proteins_contigid.txt  combined_transcripts.fasta  | seqkit seq -w 0  > 01_MetaGeneMark_Prodigal_predicted_proteins_contigid.fasta

sort all_contigid.txt  > all_contigid_sorted.txt

comm -2 -3 all_contigid_sorted.txt  01_MetaGeneMark_Prodigal_predicted_proteins_contigid.txt  >  02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.txt

seqkit grep -f 02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.txt  combined_transcripts.fasta  | seqkit seq -w 0  > 02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.fasta

samtools  faidx  01_MetaGeneMark_Prodigal_predicted_proteins_contigid.fasta

samtools  faidx  02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.fasta

cut -f 2  01_MetaGeneMark_Prodigal_predicted_proteins_contigid.fasta.fai | sort -n | awk '{print "Contig_"NR "\t" $0}'  | sed '1i Name\tLength'  >   01_MetaGeneMark_Prodigal_predicted_proteins_contigs_length.txt

cut -f 2  02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.fasta.fai  | sort -n | awk '{print "Contig_"NR "\t" $0}'   | sed '1i Name\tLength' >   02_MetaGeneMark_Prodigal_unpredicted_proteins_contigs_length.txt
```

### **uniref90_annotation**

```
sbatch uniref90_annotation.slurm
```

uniref90_annotation.slurm

```
#!/bin/bash
#SBATCH --job-name=uniref90_annotation        # job 名字
#SBATCH --output=%j.out            # 标准输出 (%j = job ID)
#SBATCH --error=%j.err             # 错误输出
#SBATCH --time=96:00:00                   # 最大运行时间，根据实际修改
#SBATCH --cpus-per-task=8                 # 使用的 CPU 核心数（与脚本里面 -t 参数一致）
#SBATCH -p Lulab


source ~/.bashrc
conda activate bio


diamond blastx \
  -d /data2/lulab1/jinyunfan/uniref/diamond/uniref90.dmnd \
  -q /lulabdata3/zhangqin/SLE_microbial_ncRNA/Prodigal_MetageneMark_combined_analysis/02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.fasta \
  -o MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_uniref90_annotation_result.txt \
  --outfmt 6 \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  --threads 9 \
  > MetaGeneMark_Prodigal_unpredicted_proteins_uniref90_annotation_diamond_run.log 2>&1
```

```
cut -f 1  MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_uniref90_annotation_result.txt  | sort | uniq  > MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_uniref90_annotation_result.id

comm -2 -3   /lulabdata3/zhangqin/SLE_microbial_ncRNA/Prodigal_MetageneMark_combined_analysis/02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.txt  MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_uniref90_annotation_result.txt.id  > MetaGeneMark_Prodigal_unpredicted_proteins_identified_not-in_uniref90_annotation_result.txt.id
```

### **nr_annotation**

``` 
sbatch nr_annotation.slurm
```

nr_annotation.slurm

```
#!/bin/bash
#SBATCH --job-name=nr_annotation        # job 名字
#SBATCH --output=%j.out            # 标准输出 (%j = job ID)
#SBATCH --error=%j.err             # 错误输出
#SBATCH --time=96:00:00                   # 最大运行时间，根据实际修改
#SBATCH --cpus-per-task=8                 # 使用的 CPU 核心数（与脚本里面 -t 参数一致）
#SBATCH -p Lulab


source ~/.bashrc
conda activate bio


diamond blastx \
  -d /data2/lulab1/zhangqin/database_nr/nr.dmnd \
  -q /lulabdata3/zhangqin/SLE_microbial_ncRNA/Prodigal_MetageneMark_combined_analysis/02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.fasta \
  -o MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_nr_annotation_result.txt \
  --outfmt 6 \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  --threads 8 \
  > MetaGeneMark_Prodigal_unpredicted_proteins_nr_annotation_diamond_run.log 2>&1
```

```
cut -f 1  MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_nr_annotation_result.txt  | sort | uniq  > MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_nr_annotation_result.id

comm -2 -3   /lulabdata3/zhangqin/SLE_microbial_ncRNA/Prodigal_MetageneMark_combined_analysis/02_MetaGeneMark_Prodigal_unpredicted_proteins_contigid.txt  MetaGeneMark_Prodigal_unpredicted_proteins_identified_in_nr_annotation_result.txt.id  > MetaGeneMark_Prodigal_unpredicted_proteins_identified_not-in_nr_annotation_result.txt.id
```

