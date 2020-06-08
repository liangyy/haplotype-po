import argparse
parser = argparse.ArgumentParser(prog='get_sample_list.py', description='''
    Create a list of sample from some input file (by extracting the columns)
''')

parser.add_argument('--input', help='''
    the file list
''')
parser.add_argument('--cols', nargs='+', help='''
    columns to use
''')
parser.add_argument('--output', help='''
    output TXT
''')
args = parser.parse_args()

import logging, time, sys
import pandas as pd

sys.path.insert(0, '../')
import misc_helper

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

df = misc_helper.read_table(args.input)
sample_list = []
for col in args.cols:
    tmp_ = df[col].tolist()
    sample_list += tmp_
with open(args.output, 'w') as f:
    for i in sample_list:
        f.write(i + '\n')

