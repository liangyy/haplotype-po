import argparse
parser = argparse.ArgumentParser(prog='generate_snp_list.py', description='''
    Generate YAML for impute_po_otf.py
''')

parser.add_argument('--yaml-template', help='''
    YAML template
''')
parser.add_argument('--trait-yaml', help='''
    trait1:
        dict_which_can_fit_into_yaml_template
''')
parser.add_argument('--out-yaml', default=None, help='''
    Path to output YAML
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

import yaml, os
os.path.insert(0, '../prs')
import gwas_reader


out_dict = {}

with open(args.yaml_template, 'r') as f:
    template_lines = f.readlines()
traits = gwas_reader.read_yaml(args.trait_yaml)

for trait in traits.keys():
    kws = traits[trait]
    trait_str = template_lines.format(**kws)
    out_dict[trait] = yaml.safe_load(trait_str)

with open(args.out_yaml, 'w') as f:
    yaml.dump(out_dict, f)

    
