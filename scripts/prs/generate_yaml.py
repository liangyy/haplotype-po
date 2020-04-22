import argparse
parser = argparse.ArgumentParser(prog='generate_yaml.py', description='''
    Generate YAML for naive_prs.py
''')

parser.add_argument('--yaml-template', help='''
    YAML template
''')
parser.add_argument('--gwas-list', help='''
    List of GWAS names to fill in: {gwas}
''')
parser.add_argument('--chr-num', default=None, help='''
    If set, will fill in chromosome number: {chr_num}
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

import yaml, copy
import gwas_reader


out_dict = {}

template = gwas_reader.read_yaml(args.yaml_template)
with open(args.gwas_list, 'r') as f:
    for gwas in f:
        gwas = gwas.strip()
        logging.info(f'Adding on {gwas}')
        out_dict[gwas] = copy.deepcopy(template)
        out_dict[gwas]['gwas']['path'] = out_dict[gwas]['gwas']['path'].format(gwas=gwas)
        if args.chr_num is None:
            out_dict[gwas]['ld_clump'] = out_dict[gwas]['ld_clump'].format(gwas=gwas)
        if args.chr_num is not None:
            out_dict[gwas]['ld_clump'] = out_dict[gwas]['ld_clump'].format(gwas=gwas, chr_num=args.chr_num)
        # print(out_dict)

with open(args.out_yaml, 'w') as f:
    yaml.dump(out_dict, f)

    
