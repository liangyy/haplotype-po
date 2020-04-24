import argparse
parser = argparse.ArgumentParser(prog='run_bgen2hdf5.py', description='''
    Convert BGEN into HDF5. 
    Work with ukb_hap_v2 BGEN.
''')

parser.add_argument('--bgen', help='''
    BGEN file path
''')
parser.add_argument('--bgi', help='''
    BGEN BGI file path
''')
parser.add_argument('--sample', help='''
    BGEN SAMPLE file path
''')
parser.add_argument('--output-hdf5', help='''
    HDF5 file name of output (if not exists, it will be created)
''')
parser.add_argument('--snp-chunk-size', type=int, default=100, help='''
    Number of SNPs to process at a time
''')
parser.add_argument('--bgen-writing-cache-size', type=int, default=50, help='''
    BGEN reading cache size in MB. (It should be set carefully: 
    pre-factor x nvariant x size of dtype x snp chunk size / 1024 ** 2
    where pre-factor ~ 10)
''')
parser.add_argument('--max-sample-chunk-size', type=int, default=10000, help='''
    Maximum size of chunk on sample axis.
''')
parser.add_argument('--max-snp-chunk-size', type=int, default=100, help='''
    Maximum size of chunk on snp axis. 
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

from tqdm import tqdm
import sys
sys.path.insert(0, '../prs')
import ukb_hap_reader
import geno_hdf5


def size_in_mb_to_cache_size(size):
    return int(size * (1024 ** 2))    

logging.info('Loading UKB HAP BGEN reader')
reader = ukb_hap_reader.UKBhapReader(
    bgen_path=args.bgen,
    bgen_bgi_path=args.bgi,
    sample_path=args.sample
)

logging.info('Building variant list')
pos = [ i.split(':')[1] for i in reader.variant_index.keys() ]
non_effect_allele = [ i.split(':')[2] for i in reader.variant_index.keys() ]
effect_allele = [ i.split(':')[3] for i in reader.variant_index.keys() ]
chrom = [ '' for i in reader.variant_index.keys() ]

logging.info('Loading sample list')
with open(sample, 'r') as f:
    sample_list = f.readlines()
    sample_list.pop(0)
    sample_list.pop(0)
sample_list = [ i.split(' ')[0] for i in sample_list ]

logging.info('Building variant generator')
generator = reader.retrieve_from_list(
    chrom=chrom,
    pos=pos,
    non_effect_allele=non_effect_allele,
    effect_allele=effect_allele,
    n_var_cached=args.snp_chunk_size
)

logging.info('Initializing genotype writer')
geno_writer = geno_hdf5.GenotypeHDF5(
    list_sample=sample_list,
    num_variants=len(pos),
    output_h5='test-bgen2hdf5.h5',
    dict_variant_meta={
        'position': 'position',
        'chromosome': 'chr',
        'reference_allele': 'allele0',
        'dosage_allele': 'allele1'
    },
    phased=True,
    dtype=int,
    cache_size=size_in_mb_to_cache_size(args.bgen_writing_cache_size)
)

logging.info('Converting')
for dosage_row in tqdm(generator, total=len(chrom)):
    geno_writer.fill_in(dosage_row, ['haplo_dosage_1', 'haplo_dosage_2'])

logging.info('Saving') 
geno_writer.save()
