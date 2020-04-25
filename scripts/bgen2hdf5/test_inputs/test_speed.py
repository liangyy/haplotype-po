# test I/O for HDF5 and BGEN
# on ANL nucleus
# ukb_hap_v2 chromosome 16

import argparse
parser = argparse.ArgumentParser(prog='test_speed.py', description='''
    Test I/O for HDF5 and BGEN
    on ANL nucleus
    with ukb_hap_v2 chromosome 16
''')
parser.add_argument('--nsnp', type=int, default=10, help='''
    Number of SNPs to load.
    Don't use number greater than 5000. 
    It will exceed the range of SNPs in test file.
''')
args = parser.parse_args()

import h5py
import sys
sys.path.insert(0, '../../../code/')
import lib_test
import imp
import time
import numpy as np
lib_test = imp.reload(lib_test)
hdf5 = '/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr16.h5'
bgen = '/vol/bmd/meliao/data/haplotype/hap/ukb_hap_chr16_v2.bgen'
bgi = '/vol/bmd/meliao/data/haplotype/hap_bgi/ukb_hap_chr16_v2.bgen.bgi'

start_pos = [ 0, 10, 100, 1000, 10000, 15000 ]
infostr = '##NSNP={}##'.format(args.nsnp)
for pos in start_pos:
    print('{} # Extracting {} -- {}'.format(infostr, pos, pos + args.nsnp))
    print(f'{infostr} # 1. Testing speed')
    # print('** HDF5')
    t0 = time.time()
    pos1, geno1 = lib_test.query_hdf5(pos, args.nsnp, hdf5)
    print(f'{infostr} # ** HDF5 takes', time.time() - t0, 'seconds', (time.time() - t0) / args.nsnp, 'per snp')
    t0 = time.time()
    var_info2, geno2 = lib_test.query(pos1.astype(int), bgen, bgi)
    print(f'{infostr} # ** RBGEN_ranges takes', time.time() - t0, 'seconds', (time.time() - t0) / args.nsnp, 'per snp')
    t0 = time.time()
    var_info3, geno3 = lib_test.query(var_info2.rsid, bgen, bgi, by='id')
    print(f'{infostr} # ** RBGEN_rsids takes', time.time() - t0, 'seconds', (time.time() - t0) / args.nsnp, 'per snp')


    print(f'{infostr} # 2. Checking if the extracted genotypes agree')
    haplo = { k: { 1: None, 2: None } for k in range(1, 4) }
    code_dic = { 1: 'HDF5', 2: 'BGEN_ranges', 3: 'RBGEN_rsids'}
    haplo[1][1] = geno1[0, :, :]  # np.einsum('ij->ji', geno1[0, :, :])
    haplo[1][2] = geno1[1, :, :]  # np.einsum('ij->ji', geno1[1, :, :])
    haplo[2][1], haplo[2][2] = lib_test.geno_prob_to_haplo(geno2)
    haplo[3][1], haplo[3][2] = lib_test.geno_prob_to_haplo(geno3)
    iffail = False
    for i in range(1, 4):
        for j in range(i + 1, 4):
            for hh in range(1, 3):
                if lib_test.check_equal(haplo[i][hh], haplo[j][hh]) is False:
                    print(f'{infostr} # Failed the checking haplotype {hh}:', code_dic[i], 'vs', code_dic[j])
                    iffail = True
                    break
    if iffail is False:
        print(f'{infostr} # ** Checking passed!')
    else:
        print(f'{infostr} # ** Checking FAILED!')

