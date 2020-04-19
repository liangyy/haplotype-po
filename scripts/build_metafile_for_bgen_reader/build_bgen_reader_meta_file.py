import argparse
parser = argparse.ArgumentParser(prog='build_bgen_reader_meta_file.py', description='''
    Build meta file for bgen_reader
''')
parser.add_argument('--bgen-pattern', help='''
    Input BGEN, for example:
    prefix.chr{chr}.suffix
''')
parser.add_argument('--out-prefix', help='''
    Output prefix
''')
args = parser.parse_args()

import logging, time, sys, os
# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

dirname = os.path.dirname(os.path.abspath(args.out_prefix))
if not os.path.isdir(dirname):
    os.mkdir(dirname)

import bgen_reader

for i in range(1, 23):
    i = str(i)
    metafile_filepath = '{outprefix}.chr{i}.metafile'.format(
        outprefix=args.out_prefix,
        i=i
    )
    logging.info(f'Processing chr{i}')
    bgen_reader.create_metafile(
        args.bgen_pattern.format(chr=i),
        metafile_filepath
    )
    

