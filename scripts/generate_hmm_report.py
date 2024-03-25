import pandas as pd
import numpy as np
import argparse
import os
from re import search


description = """Script make reports file from a result file of script 'analyze_seq.py'"""


parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i', '--input', type=str, help='File or regular expression of file')
parser.add_argument('-i2', '--input2', type=str, help='File or regular expression of file')
parser.add_argument('-o', '--output', type=str, help='Output folder')
parser.add_argument('-m', '--hmm_info', type=str, help='A folder with drawings')
parser.add_argument('-t', '--hmm_out', type=str, help='A one column csv with models, that have passed the treshold')
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
    
if args.hmm_out:
    hmm_output = args.hmm_out
else:
    hmm_output = False

if args.hmm_info:
    if_condition(os.path.isfile(args.hmm_info), "drawings folder doesn't exist")
    hmm_info = args.hmm_info
    hmm_info = pd.read_csv(hmm_info, index_col=0)
else:
    if_condition(False, 'Missed -h argument')

    
def export_seq_data(seq):
    name, chain, part, length, length_contig, skew, sdf = seq.split(';')
    chain, part, length, skew, length_contig = chain.replace('chain_', ''), part.replace('part_', ''), length.replace('len_', ''), skew.replace('skew_', ''), length_contig.replace('len_contig_', '')
    frame = int(chain.replace('+1', '1').replace('-1', '4')) + int(skew)
    return [name, part, length, frame, length_contig]


def get_coordinate(dataframe):
    result = []
    for i, row in dataframe.iterrows():
        if row['Frame'] <= 3:
            from_ = row['From'] * 3 + row['Part']
            to_ = row['To'] * 3 + row['Part']
            result.append([from_, to_])
        elif row['Frame'] >= 4:
            from_ = (row['Lengh'] - row['From'] * 3)  + row['Part']
            to_ = (row['Lengh'] - row['To'] * 3)  + row['Part']
            result.append([from_, to_])
    return np.array(result)
                    


def create_report(input_file, hmm_info):
    report = input_file.merge(hmm_info, left_on='Query', right_on='Model_ID', how='left').drop('Model_ID', axis=1)
    report['Score_ratio'] = report.Score/report.Threshold
    report = report.sort_values('Score_ratio', ascending=False)
    try:
        report[['Name', 'Part', 'Lengh', 'Frame', 'Length_contig']] = np.array(report.Name.apply(lambda  x: export_seq_data(x)).to_list())
        report[['Lengh', 'From', 'To', 'Part', 'Frame', 'Length_contig']] = report[['Lengh', 'From', 'To', 'Part', 'Frame', 'Length_contig']].astype('int')
        report[['From', 'To']] = get_coordinate(report)
        report = report.drop('Part', axis=1)
    except Exception:
        None
    return report


report = create_report(input_file, hmm_info)

if args.input2:
    report2 = create_report(input_file2, hmm_info)
    report = pd.concat([report, report2])
    


    
report = report.query('Score_ratio > 1.5')

if hmm_output:
    report[['Query']].drop_duplicates().to_csv(hmm_output, header=None, index=None)

report.reset_index(drop=True).to_csv(output)

print(f'Successfully end with file "{args.input}"')
