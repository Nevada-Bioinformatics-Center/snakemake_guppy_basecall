import glob

configfile: "config.yaml"

inputdirectory=config["directory"]
outdirectory=config["directory"]
#BARCODES, SAMPLES, = glob_wildcards(inputdirectory+"/fast5_pass/{barcode}/{sample}.fast5", followlinks=True)
#print("Pass Sample List")
#print(SAMPLES)
#print(BARCODES)
#
#SAMPLES_skip, = glob_wildcards(inputdirectory+"/fast5_skip/{sample_skip}.fast5", followlinks=True)
#print("Skip Sample List")
#print(SAMPLES_skip)

SAMPLES_gen, = glob_wildcards(inputdirectory+"/fast5/{sample_gen}.fast5", followlinks=True)
print("General Sample List")
print(SAMPLES_gen)

#wildcard_constraints:
#    barcode="barcode\d+"
#    sample="\w+\d+_\w+_\w+\d+_.+_\d"


##### target rules #####
rule all:
    input: 
#       expand(inputdirectory+'/basecall/pass/{sample}_{barcode}/', zip, barcode=BARCODES, sample=SAMPLES),
#       expand(inputdirectory+'/basecall/skip/{sample_skip}/', sample_skip=SAMPLES_skip),
       expand(inputdirectory+'/basecall/output/{sample_gen}/', sample_gen=SAMPLES_gen),
       directory=directory(outdirectory+"/basecall/demultiplex/"),


#rule make_indvidual_samplefiles_pass:
#    input:
#        inputdirectory+"/fast5_pass/{barcode}/{sample}.fast5",
#    output:
#        "lists/{sample}_{barcode}.txt",
#    shell:
#        "basename {input}  > {output}"
#
#rule make_indvidual_samplefiles_skip:
#    input:
#        inputdirectory+"/fast5_skip/{sample_skip}.fast5",
#    output:
#        "lists/{sample_skip}.txt",
#    shell:
#        "basename {input}  > {output}"

rule make_indvidual_samplefiles_gen:
    input:
        inputdirectory+"/fast5/{sample_gen}.fast5",
    output:
        "lists/{sample_gen}.txt",
    shell:
        "basename {input}  > {output}"


#rule guppy_basecall_persample_pass:
#    input:
#        directory=inputdirectory+"/fast5_pass/{barcode}/",
#        samplelist="lists/{sample}_{barcode}.txt",
#    output:
#        directory=directory(outdirectory+"/basecall/pass/{sample}_{barcode}/"),
#    params: 
#        basealgo=config["basealgo"],
#    shell:
#        "guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params.basealgo} -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"
#
#rule guppy_basecall_persample_skip:
#    input:
#        directory=inputdirectory+"/fast5_skip/",
#        samplelist="lists/{sample_skip}.txt",
#    output:
#        directory=directory(outdirectory+"/basecall/skip/{sample_skip}/"),
#    params: 
#        basealgo=config["basealgo"],
#    shell:
#        "guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params.basealgo} -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"
#        #"guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params.basealgo} --compress_fastq -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200"

rule guppy_basecall_persample:
    input:
        directory=inputdirectory+"/fast5/",
        samplelist="lists/{sample_gen}.txt",
    output:
        directory=directory(outdirectory+"/basecall/output/{sample_gen}/"),
    params: 
        #basealgo=config["basealgo"],
        kit=config["kit"],
        flowcell=config["flowcell"],
    shell:
        "guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} --flowcell {params.flowcell} --kit {params.kit} -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200 --chunk_size 2000"
        #"guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} -c {params.basealgo} -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200 --chunk_size 500"
        #"guppy_basecaller -i {input.directory} --input_file_list {input.samplelist} -s {output.directory} --flowcell {params.flowcell} --kit {params.kit} --barcode_kits \"{params.kit}\" --trim_barcodes --trim_adapters --trim_primers -x \"auto\" --gpu_runners_per_device 3 --num_callers 2 --chunks_per_runner 200 --chunk_size 500"

rule guppy_barcoder:
    input:
        expand(outdirectory+"/basecall/output/{sample_gen}/", sample_gen=SAMPLES_gen),
        indir=outdirectory+"/basecall/output/",
    output:
        directory=directory(outdirectory+"/basecall/demultiplex/"),
    params: 
        kit=config["kit"],
    shell:
        "guppy_barcoder -i {input.indir} -s {output.directory} --barcode_kits {params.kit} -x \"auto\" --trim_barcodes --trim_adapters --trim_primers --recursive --fastq_out --records_per_fastq 0"
