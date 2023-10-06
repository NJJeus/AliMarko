import pandas as pd
import numpy as np
import argparse
import os


description = """Script make reports file from a result file of script 'analyze_seq.py'"""


parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i', '--input', type=str, help='File or regular expression of file')
parser.add_argument('-i2', '--input2', type=str, help='File or regular expression of file')
parser.add_argument('-o', '--output', type=str, help='Output folder')
parser.add_argument('-m', '--hmm_info', type=str, help='A folder with drawings')
args = parser.parse_args()

def if_condition(x, message):
    if not x:
        print(message)
        exit()


if args.input:
    if_condition(os.path.isfile(args.input), "An input file does not exists")
    input_file = args.input
    input_file = pd.read_csv(input_file, index_col=0)
else:
    if_condition(False, 'Missed -i argument')
    
if args.input2:
    if_condition(os.path.isfile(args.input2), "An input file 2 does not exists")
    input_file2 = args.input2
    input_file2 = pd.read_csv(input_file2, index_col=0)


if args.output:
    output = args.output
else:
    if_condition(False, 'Missed -o argument')

if args.hmm_info:
    if_condition(os.path.isfile(args.hmm_info), "drawings folder doesn't exist")
    hmm_info = args.hmm_info
    hmm_info = pd.read_csv(hmm_info, index_col=0)
else:
    if_condition(False, 'Missed -h argument')




def create_report(input_file, hmm_info):
    report = input_file.merge(hmm_info, left_on='Query', right_on='Model_ID', how='left').drop('Model_ID', axis=1)
    report['Score_ratio'] = report.Score/report.Threshold
    report = report.sort_values('Score_ratio', ascending=False)
    
    return report


report = create_report(input_file, hmm_info)

if args.input2:
    report2 = create_report(input_file2, hmm_info)
    report = pd.concat([report, report2])

    
report = report.query('Score_ratio > 1.5')#.drop('Name', axis=1)
#report = report.groupby(['Query', 'Taxon']).mean().reset_index().sort_values('Score_ratio', ascending=False)

report.reset_index(drop=True).to_csv(output)

print(f'Successfully end with file "{args.input}"')