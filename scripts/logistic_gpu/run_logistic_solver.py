import argparse
parser = argparse.ArgumentParser(prog='run_logistic_solver.py', description='''
    Run logistic regression on a set of variants.
''')

parser.add_argument('', help='''
''')
args = parser.parse_args()

import logging, time, sys
# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

TODO
