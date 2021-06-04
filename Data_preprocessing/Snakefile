import os
import fnmatch

########################################################
# global varibales
# 1. modify global varibales here
# 2. modify cluster.json
########################################################

rawData_dir = "./raw_data" # raw sequencing data directory
processedData_dir = './Processed_data' # processed sequencing data directory
fastq_ext = '.fq.gz' # fastq file extension

core = 20 # threads to use
core_cellRanger = int(core*0.8) # threads to use in cellRanger
mem = 64 # memory to use in cellRanger

bbduk_dir = "./bbmap/bbduk.sh" # bbduk.sh directory
cellranger_dir = "./cellranger-atac-1.2.0/cellranger-atac" # cellranger-atac directory
ref_dir = "./refdata-cellranger-atac-mm10-1.2.0" # reference data directory

########################################################

########################################################
# Do not change codes below
########################################################

samples = os.listdir(rawData_dir)
output = [processedData_dir+'/'+sample+'/CellRanger/'+sample+'/outs/fragments.tsv.gz' for sample in samples]

list_fastq={}
for sample in samples:
    list_fastq[sample]=(fnmatch.filter(os.listdir(rawData_dir+'/'+sample), '*1'+fastq_ext)[0].split(".")[0],
                        fnmatch.filter(os.listdir(rawData_dir+'/'+sample), '*2'+fastq_ext)[0].split(".")[0])
                        
# generating folders

for sample in samples:
    
    output_dir=processedData_dir+'/'+sample+'/CellRanger'
    tmp_data=output_dir+'/tmp_data'
    qc_raw_data=tmp_data+'/qc_raw_data'
    log_dir=output_dir+'/log'
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(tmp_data):
        os.makedirs(tmp_data)
    if not os.path.exists(qc_raw_data):
        os.makedirs(qc_raw_data)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)


rule all:
  input:
    output

rule filter_L1:
  input:
    in1 = lambda wildcards: rawData_dir+"/{sample}/"+list_fastq[wildcards.sample][0]+fastq_ext,
    in2 = lambda wildcards: rawData_dir+"/{sample}/"+list_fastq[wildcards.sample][1]+fastq_ext
  output:
    out1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R1.fastq.gz",
    out2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R2.fastq.gz"
  shell:
    '''
    {bbduk_dir} \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k=34 mm=f rcomp=f restrictleft=85 skipr2=t \
    hdist=3 \
    stats={processedData_dir}/{wildcards.sample}/CellRanger/tmp_data/qc_raw_data/{wildcards.sample}_stats.linker1.txt \
    threads={core} \
    literal=AGATGTGTATAAGAGACAGCATCGGCGTACGACT
    '''
    
rule filter_L2:
  input:
    in1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R1.fastq.gz",
    in2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R2.fastq.gz"
  output:
    out1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R1.fastq.gz",
    out2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R2.fastq.gz"
  shell:
    '''
    {bbduk_dir} \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k=30 mm=f rcomp=f restrictleft=43 skipr2=t \
    hdist=3 \
    stats={processedData_dir}/{wildcards.sample}/CellRanger/tmp_data/qc_raw_data/{wildcards.sample}_stats.linker2.txt \
    threads={core} \
    literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG
    '''

rule bc_process:
  input:
    processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R1.fastq.gz"
  output:
    out1 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R1_001.fastq",
    out2 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R2_001.fastq"
  shell:
    '''
    python BC_process.py --input {input} --output_R1 {output.out1} --output_R2 {output.out2}
    '''
    
rule R2_rename:
  input:
    processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R2.fastq.gz"
  output:
    processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R3_001.fastq.gz"
  shell:
    '''
    cp {input} {output}
    '''
    
rule cell_ranger:
  input:
    in1 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R1_001.fastq",
    in2 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R2_001.fastq",
    in3 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R3_001.fastq.gz"
  output:
    processedData_dir + "/{sample}/CellRanger/{sample}/outs/fragments.tsv.gz"
  shell:
    '''
    cd {processedData_dir}/{wildcards.sample}/CellRanger/
    rm -r {wildcards.sample}
    
    {cellranger_dir} count \
    --id={wildcards.sample} \
    --reference={ref_dir} \
    --fastqs={processedData_dir}/{wildcards.sample}/CellRanger/tmp_data \
    --sample={wildcards.sample} \
    --localcores={core_cellRanger} \
    --localmem={mem}
    '''
    