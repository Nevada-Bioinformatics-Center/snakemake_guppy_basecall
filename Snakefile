import glob

configfile: "config.yaml"

inputdirectory=config["directory"]
outdirectory=config["directory"]
BARCODES, SAMPLES, = glob_wildcards(inputdirectory+"/fast5_pass/{barcode}/{sample}.fast5", followlinks=True)
print("Pass Sample List")
print(SAMPLES)
print(BARCODES)

SAMPLES_skip, = glob_wildcards(inputdirectory+"/fast5_skip/{sample_skip}.fast5", followlinks=True)
print("Skip Sample List")
print(SAMPLES_skip)

#wildcard_constraints:
#    barcode="barcode\d+"
#    sample="\w+\d+_\w+_\w+\d+_.+_\d"


##### target rules #####
rule all:
    input: 
       expand(inputdirectory+'/basecall/pass/{sample}_{barcode}/', zip, barcode=BARCODES, sample=SAMPLES),
       expand(inputdirectory+'/basecall/skip/{sample_skip}/', sample_skip=SAMPLES_skip),


rule make_indvidual_samplefiles_pass:
    input:
        inputdirectory+"/fast5_pass/{barcode}/{sample}.fast5",
    output:
        "lists/{sample}_{barcode}.txt",
    shell:
        "basename {input}  > {output}"

rule make_indvidual_samplefiles_skip:
    input:
        inputdirectory+"/fast5_skip/{sample_skip}.fast5",
    output:
        "lists/{sample_skip}.txt",
    shell:
        "basename {input}  > {output}"


rule guppy_basecall_persample_pass:
    input:
        directory=inputdirectory+"/fast5_pass/{barcode}/",
        samplelist="lists/{sample}_{barcode}.txt",
    output:
        directory=directory(outdirectory+"/basecall/pass/{sample}_{barcode}/"),
    params: 
        basealgo=config["basealgo"],
    shell:
        "guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params.basealgo} -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"

rule guppy_basecall_persample_skip:
    input:
        directory=inputdirectory+"/fast5_skip/",
        samplelist="lists/{sample_skip}.txt",
    output:
        directory=directory(outdirectory+"/basecall/skip/{sample_skip}/"),
    params: 
        basealgo=config["basealgo"],
    shell:
        "guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params.basealgo} -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"
        #"guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params.basealgo} --compress_fastq -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"


#rule guppy_linkfastq:
#    input:
#        #glob_wildcards("basecall/{sample}/*/*.fastq.gz"),
#        "basecall/{sample}/pass/*.fastq.gz",
#    output:
#        "basecall/{sample}.fastq.gz",
#    shell:
#        "ln -s {input} {output}"
#
#rule fastqc_pretrim:
#    input:
#        #"basecall/{sample}/{failpass}/{runid}.fastq.gz",
#        "basecall/{sample}.fastq.gz"
#    output:
#        html="qc/fastqc_pretrim/{sample}.html",
#        zip="qc/fastqc_pretrim/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#    params: ""
#    log:
#        "logs/fastqc_pretrim/{sample}.log"
#    threads: 1
#    wrapper:
#        "v0.75.0/bio/fastqc"
#
#rule multiqc:
#    input:
#        #expand("basecall/{sample}.fastq.gz", sample=SAMPLES)
#        expand("qc/fastqc_pretrim/{sample}_fastqc.zip", sample=SAMPLES)
#    output:
#        "qc/multiqc.html"
#    params:
#        ""  # Optional: extra parameters for multiqc.
#    log:
#        "logs/multiqc.log"
#    wrapper:
#        "0.77.0/bio/multiqc"

#rule fastqc_pretrim:
#    input:
#        "basecall/{sample}/{failpass}/{runid}.fastq.gz",
#    output:
#        html="qc/fastqc_pretrim/{sample}_{failpass}_{runid}.html",
#        zip="qc/fastqc_pretrim/{sample}_{failpass}_{runid}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#    params: ""
#    log:
#        "logs/fastqc_pretrim/{sample}_{failpass}_{runid}.log"
#    #resources: time_min=320, mem_mb=8000, cpus=1
#    threads: 1
#    wrapper:
#        "v0.75.0/bio/fastqc"
