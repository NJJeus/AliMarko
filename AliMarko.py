import subprocess
import argparse
import os

description = """AliMarko wrapper"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-m', '--paired', type=str, help='')
parser.add_argument('-s', '--suffix', type=str, help='')
parser.add_argument('-i', '--input_directory', type=str, help='', default="DATA_test/raw_fastq/")
parser.add_argument('-o', '--output_directory', type=str, help='', default="DATA_test4/")
parser.add_argument('-g', '--genome_reference', type=str, help='', default='ictv_virus_reference.fa')
parser.add_argument('-f', '--hmm_folder', type=str, help='', default='test_HMM/')
parser.add_argument('-n', '--hmm_info', type=str, help='', default='ictv_tables/hmm_info.csv')
parser.add_argument('-k', '--kraken_database', type=str, help='', default='16S_Greengenes_k2db')
parser.add_argument('-t', '--threads', type=int, help='', default=4)

args = parser.parse_args()

conf = vars(args)

paired = conf.pop('paired')
threads = conf.pop('threads')
current_dir = os.getcwd()

if paired == "paired":
    cmd = f'snakemake -s scripts/AliMarko_paired_snakemake.py --stats stats.json --use-conda -j {threads}  --printshellcmds --config '
elif paired == "single":
    cmd = f'snakemake -s scripts/AliMarko_single_snakemake.py --use-conda -j {threads} --printshellcmds --config --printshellcmds '
else:
    print('LOL?')
    exit()

cmd += f'workdir="{current_dir}" '
    
for arg, value in conf.items():
    cmd = cmd + f'{arg}="{value}" '

result = subprocess.run(cmd, shell=True, check=True, text=True)
