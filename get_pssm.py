import os
import re
import argparse
import codecs  # 或者io，使用哪种包无所谓
import pandas as pandas


def PSSM(pdb_name, chain, ipssm, out_pssm, blast):
    # with open('./one_pdb/fasta/' + pdb_name + '_fasta.txt','')

    os.system(
        r'psiblast -query "' + ipssm
        + str(pdb_name) + '-' + str(chain)
        + '.fa"'
        + r' -db ' + blast + 'swissprot -evalue 0.001 -num_iterations 3'
        + r' -out_ascii_pssm "' + out_pssm +
        str(pdb_name) + '-' + str(chain) + '.pssm"')

    print(r'psiblast -query ' + ipssm
          + str(pdb_name) + '-' + str(chain)
          + '.fa"'
          + r' -db ' + blast + 'swissprot -evalue 0.001 -num_iterations 3'
          + r' -out_ascii_pssm ' + out_pssm +
          str(pdb_name) + '-' + str(chain) + '.pssm"')

    fw = open(out_pssm + 'no_pssm_02.txt', 'a')

    a = os.path.exists(out_pssm + pdb_name + '-' + chain + '.pssm')
    print(out_pssm + pdb_name + '-' + chain + '.pssm')
    print(a)
    if a == False:
        with open(ipssm + str(pdb_name) + '-' + str(chain) + '.fa' )as f:
            line = f.readlines()
        fw.writelines(line)

    print(str(pdb_name) + " PSSM over")

def main():
    parser = argparse.ArgumentParser(
        description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-input', dest='input', type=str,
                        help='Specify the location of file that contain the one_pdb chains i.e. data_list.txt.txt for the specific ion of interest',
                        required=True)
    parser.add_argument('-blast', dest='blast', type=str, help='Specify the path to blast', required=True)
    parser.add_argument('-in-pssm', dest='ipssm', type=str, help='Specify the path to pssm input', required=True)
    parser.add_argument('-out-pssm', dest='out_pssm', type=str, help='Specify the path to output for pssm',
                        required=True)

    args = parser.parse_args()

    input = args.input
    ipssm = args.ipssm
    blast = args.blast
    out_pssm = args.out_pssm

    os.makedirs(ipssm, exist_ok=True)
    os.makedirs(out_pssm, exist_ok=True)

    with open(input, 'r') as fp:
        for line in fp:
            row = line.strip().split('_')
            pdb_name = row[0]
            chain = row[1]
            PSSM(pdb_name, chain, ipssm, out_pssm, blast)

    print("Done extracting PSSM info")


# PSSM('6WPS')
if __name__ == '__main__':
    main()
