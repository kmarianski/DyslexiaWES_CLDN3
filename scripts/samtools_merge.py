import os
import subprocess

threads = 32

source = '/results/bwa/'
destination = '/results/merged_bams/'

to_be_merged = {}
for file in os.listdir(source):
    sample = file.split('_')[0]
    if sample not in to_be_merged:
        to_be_merged[sample] = []
    to_be_merged[sample].append(os.path.join(source, file))

def merge(input, destination, threads=16):
    for key, val in input.items():
        if len(val) != 1:
            bashCommand = f"samtools merge -@ {threads} -o {destination}{key}.bam {' '.join(val)}"   
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        else:
            if not os.path.exists(destination+key+'.bam'):
                os.symlink(val[0], destination+key+'.bam')

merge(to_be_merged, destination, threads)