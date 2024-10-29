# fw = open('./data/fasta/570_fasta_one_model_total.fa', 'w')
#此版本，两列有一个 是 FE2 就可以

import argparse
import os

parser = argparse.ArgumentParser(
    description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
# parser.add_argument('-res', dest='res', type=str, help='Specify the path to the  surface', required=True)
parser.add_argument('-ion', dest='ion', type=str, help='Specify the path to the features', required=True)
parser.add_argument('-input', dest='datalist_path', type=str, help='Specify the path to the features', required=True)
parser.add_argument('-label-path', dest='label_path', type=str, help='Specify the path to the features', required=True)
parser.add_argument('-pos-file-path', dest='pos_file_path', type=str, help='Specify the path to the features', required=True)
parser.add_argument('-pdb-path', dest='pdb_path', type=str, help='Specify the path to the features', required=True)

args = parser.parse_args()

# res = args.res
ion = args.ion
label_path = args.label_path
pos_file_path = args.pos_file_path
os.makedirs(label_path, exist_ok=True)

datalist_path = args.datalist_path
pdb_path = args.pdb_path

candidate_residues={
                    'ZN': ['C', 'H', 'E', 'D'],
                    'CU': ['C', 'H'],
                    'FE2': ['D', 'E', 'H'],
                    'CA': ['D', 'E', 'G', 'N'] ,
                    'MG': ['D', 'E', 'N'],
                    'MN': ['D','E','H'],
                    'NA': ['D','G','N','S'],
                    'K': ['D','E','G','N','S'],
                   }

fw = open(pos_file_path, 'w')

res=['GLY','ALA','VAL','LEU','ILE','PRO','PHE','TYR','TRP','SER','THR','CYS','MET','ASN','GLN','ASP','GLU','LYS','ARG','HIS']

with open(datalist_path, 'r') as fp:
    for lines in fp:
        row = lines.strip().split()
        pdb = row[0][:4].upper()
        chain = row[0][5]
        # print('new_postive: ' + pdb)
        # fw.write('>'+one_pdb.upper() + '-' + chain + '\n')

        # ff = open('./data/fasta/' + str(one_pdb) + '-' + str(chain) + '.fa','w')
        # ff.write('>'+one_pdb.upper() + '-' + chain + '\n')
        ion = ion
        residues = candidate_residues[ion]
        with open(f'{pdb_path}{pdb}.pdb', 'r') as f:
            for line1 in f:
                row = line1.strip().split()
                if (row[0] == 'LINK'):
                    # if (line1.find(ion) != -1):
                    if (len(row[1]) == 7):
                        # row[9] = row[8]
                        row[8] = row[7]
                        row[7] = row[6]
                        row[6] = row[5]
                        row[5] = row[4]
                        # row[5] = row[4]
                        row[4] = row[3]
                        row[3] = row[2]
                        row[2] = row[1][-3:]
                        # row[1] = row[1][:3]
                    if (len(row[2]) == 1 and row[2].isdigit() == False and row[2] != ion):  # LINK FE A FE or LINK K K
                        row.append('0')
                        row[2] = row[3]
                        row[3] = row[4]
                        row[4] = row[5]
                        row[5] = row[6]
                        row[6] = row[7]
                        row[7] = row[8]
                        row[8] = row[9]
                        row[9] = row[10]
                    if (len(row[3]) != 1):
                        row[8] = row[7]
                        row[7] = row[6]
                        row[6] = row[5]
                        row[5] = row[4]
                            # row[5] = row[4]
                        row[4] = row[3][1:]
                        row[3] = row[3][0]
                    if (len(row[5]) == 7):
                            # row.append('0')
                            # row[11] = row[10]
                        row[8] = row[7]
                        row[7] = row[6]
                        row[6] = row[5][-3:]
                            # row[6] = row[7][:-3]
                    if (len(row[6]) == 1 and row[6].isdigit() == False and row[6] != ion):  # LINK FE A FE or LINK K K
                        row.append('0')
                        row[6] = row[7]
                        row[7] = row[8]
                        row[8] = row[9]
                        row[9] = row[10]
                    if (len(row[7]) != 1):
                        row[8] = row[7][1:]
                        row[7] = row[7][0]
                    if (len(row[2]) != 1):
                        row[2] = row[2][-3:]
                    if (len(row[6]) != 1):
                        row[6] = row[6][-3:]

                        # 两列值有一列是FE2 算 label -------------------------------------------------------------
                        # if row[5] == ion and row[5] == row[6]:
                    if row[5] == ion or row[6] == ion:
                        m = row[2]  # 氨基酸
                        c = row[3]  # 链
                        c_i = 3
                    else:
                            # if row[1] == ion and row[1] == row[2]:  #
                        if row[1] == ion or row[2] == ion:
                            # if row[1] in res or row[2] in res:
                            m = row[6]
                            c = row[7]
                            c_i = 7
                        else:
                            continue  # 第一列不是ion，则结束
                    if c == chain:
                        if m == 'GLY':
                            m = 'G'
                        elif m == 'ALA':
                            m = 'A'
                        elif m == 'VAL':
                             m = 'V'
                        elif m == 'LEU':
                            m = 'L'
                        elif m == 'ILE':
                            m = 'I'
                        elif m == 'PRO':
                            m = 'P'
                        elif m == 'PHE':
                            m = 'F'
                        elif m == 'TYR':
                            m = 'Y'
                        elif m == 'TRP':
                            m = 'W'
                        elif m == 'SER':
                             m = 'S'
                        elif m == 'THR':
                                m = 'T'
                        elif m == 'CYS':
                            m = 'C'
                        elif m == 'MET':
                            m = 'M'
                        elif m == 'ASN':
                            m = 'N'
                        elif m == 'GLN':
                            m = 'Q'
                        elif m == 'ASP':
                            m = 'D'
                        elif m == 'GLU':
                            m = 'E'
                        elif m == 'LYS':
                            m = 'K'
                        elif m == 'ARG':
                            m = 'R'
                        elif m == 'HIS':
                            m = 'H'
                        if m in residues:
                            pos = row[c_i + 1]
                            fw.write(str(pdb) + '_' + str(chain) + ' ' + m + ' ' + str(c) + ' ' + str(pos) + '\n')

    # fw.write('\n')
    # ff.write('\n')

# --------------------------------------------------------------------biopython

# from Bio import SeqIO
# pdbfile = './one_pdb/antigen/7K43-E.one_pdb'
# with open(pdbfile) as handle:
#     sequence = next(SeqIO.parse(handle, "one_pdb-atom"))
# with open("one_pdb/fasta/7K43_biopython.fa", "w") as output_handle:
#     SeqIO.write(sequence, output_handle, "fasta")
# import os
#
# fa_strr=''
# fasta_path = './one_pdb/fasta/'
# fw = open('./one_pdb/fasta/7K43_biopython.fasta','a')
# with open(os.path.join(fasta_path,'7K43_biopython.fa'),'r')as f:
#     strr=f.read()
#     print(strr[0:7])
#     for i in strr[8:]:
#         if i !='\n' and i !='X':
#             fa_strr=fa_strr+i
# # print(fa_strr)
# fw.write('>7K43-E' + '\n')
# fw.write(fa_strr + '\n')
# # fw.write('\n')
