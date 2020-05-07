import argparse
parser = argparse.ArgumentParser(prog='add_baseline.py', description='''
    Add baseline to imputed prob z table.
''')

parser.add_argument('--input', help='''
    TSV.GZ input
''')
parser.add_argument('--output', help='''
    TSV.GZ output
''')
args = parser.parse_args()
import pandas as pd

df = pd.read_csv(args.input, sep='\t', compression='gzip')
df['baseline'] = 0.5
df['flip'] = 1 - df['prob_z']
df.to_csv(args.output, sep='\t', compression='gzip', index=0)

