import argparse
import os

parser = argparse.ArgumentParser(
    description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
parser.add_argument('-ipath', dest='ipath', type=str, help='Specify the path to  the specific ion of interest',
                    required=True)
parser.add_argument('-label', dest='label', type=str, help='Specify the directory for the ion under consideration',
                    required=True)
parser.add_argument('-surface', dest='surface', type=str, help='Specify the path to the  surface', required=True)
parser.add_argument('-feature', dest='feature', type=str, help='Specify the path to the features', required=True)

args = parser.parse_args()

ipath = args.ipath
label = args.label
surface = args.surface
feature = args.feature

os.makedirs(feature, exist_ok=True)


fw = open('./data/ZN/negivate_subset_2024.txt', 'w')

res=['GLY','ALA','VAL','LEU','ILE','PRO','PHE','TYR','TRP','SER','THR','CYS','MET','ASN','GLN','ASP','GLU','LYS','ARG','HIS']

with open('./data/ZN/data_list_final.txt', 'r') as fp:
    for lines in fp:
        row = lines.strip().split()
        pdb = row[0][:4].upper()
        chain = row[0][5]
        # fw.write('>'+one_pdb.upper() + '-' + chain + '\n')

        # ff = open('./data/fasta/' + str(one_pdb) + '-' + str(chain) + '.fa','w')
        # ff.write('>'+one_pdb.upper() + '-' + chain + '\n')
        ion = 'ZN'
        residues = ['C', 'H', 'E', 'D']
        with open('./data/ZN/pdb_update/{}.pdb'.format(pdb, chain), 'r') as f:
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
                        if (len(row[2]) == 1 and row[2].isdigit() == False and row[
                            2] != ion):  # LINK FE A FE or LINK K K
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
                        if (len(row[6]) == 1 and row[6].isdigit() == False and row[
                            6] != ion):  # LINK FE A FE or LINK K K
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
                        # if row[5] == ion or row[6] == ion:
                        if row[5] in res or row[6] in res:
                            m = row[2]  # 氨基酸
                            c = row[3]  # 链
                            c_i = 3
                        else:
                            # if row[1] == ion and row[1] == row[2]:  #
                            # if row[1] == ion or row[2] == ion:
                            if row[1] in res or row[2] in res:
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
