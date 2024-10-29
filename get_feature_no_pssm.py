"""
   Generate feature files (without pssm)
   python 3.7

   1. input
      all label (new_label.txt)
      all surface (new_surface1.txt)

   2. output
      test_data_feature_surface.txt
"""

import numpy as np
from rdkit import Chem
import os
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-ipath', dest='ipath', type=str, help='Specify the path to  the specific ion of interest',
                        required=True)
    parser.add_argument('-label', dest='label', type=str, help='Specify the directory for the ion under consideration',
                        required=True)
    parser.add_argument('-surface', dest='surface', type=str, help='Specify the path to the  surface', required=True)
    parser.add_argument('-feature', dest='feature', type=str, help='Specify the path to the features', required=True)
    parser.add_argument('-one-pdb-path', dest='one_pdb_path', type=str, help='Specify the directory for the ion under consideration', required=True)

    args = parser.parse_args()

    ipath = args.ipath
    label = args.label
    surface = args.surface
    feature = args.feature
    one_pdb_path = args.one_pdb_path

    os.makedirs(feature, exist_ok=True)
    # os.makedirs(feature, exist_ok=True)

    ff = open(feature + 'feature_no_pssm.txt', 'w')

    label_file = label + 'label.txt'

    surface_file = surface + 'surface.txt'

    def onek_encoding_unk(x, allowable_set):
        if x not in allowable_set:
            x = allowable_set[-1]
        return [x == s for s in allowable_set]

    with open(label_file, 'r') as num_fp:
        data_num = len(num_fp.readlines())
    point = 1
    with open(label_file, 'r') as fp, open(surface_file) as fp3:
        ff.write(str(int(data_num / 2)))
        ff.write('\t')
        ff.write('\n')
        for line1 in fp:  # line1:str '>3jcx-A\n'
            # if line1[1:5] == '1DZB':
            #     print()
            if line1.startswith('>'):
                print(line1)

                fp3.readline()
                row1 = line1.strip().split()  # list:['>3jcx-A']
                file_name = row1[0][1:5]  # 3jcx
                # file_name = file_name.lower()
                seq = row1[0][6]
                file = one_pdb_path + file_name + '_' + seq + '.pdb'
                print("get_feature_no_pssm test: " + file)
                # file = 'pdb/' + file_name + '.pdb'
                try:
                    if (Chem.MolFromPDBFile(file, removeHs=False, flavor=1, sanitize=False)):

                        fasta = fp.readline().strip()  # label
                        surface = fp3.readline().strip()  # surface
                        seq = row1[0][6]  # chain

                        mol = Chem.MolFromPDBFile(file, removeHs=False, flavor=1, sanitize=False)

                        natoms = mol.GetNumAtoms()  # atoms total in a chain

                        atom_feature = []

                        # matrix 邻接矩阵
                        # matrix = np.zeros((natoms, natoms),dtype='float32')

                        ELEM_LIST = ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'Al', 'I', 'B',
                                     'K', 'Se', 'Zn', 'H', 'Cu', 'Mn', 'unknown']

                        atom_info = []
                        i = 0
                        n = -100
                        p = -1
                        first = 0
                        last = 0
                        with open(file, 'r') as f:  # pdb
                            m = 0
                            for line in f:
                                m += 1
                                row = line.strip().split()
                                if (row[0] == 'ATOM'):
                                    if (len(row[2]) == 7):
                                        row[7] = row[6]
                                        row[6] = row[5]
                                        row[5] = row[4]
                                        row[4] = row[3]
                                        row[3] = row[2][4:]
                                        row[2] = row[2][:4]
                                    if (len(row[2]) == 8):
                                        row[7] = row[6]
                                        row[6] = row[5]
                                        row[5] = row[4]
                                        row[4] = row[3]
                                        row[3] = row[2][5:]
                                        row[2] = row[2][:5]
                                    if (len(row[4]) != 1):
                                        # row.append('0')
                                        row[8] = row[7]
                                        row[7] = row[6]
                                        row[6] = row[5]
                                        row[5] = row[4][1:]
                                        row[4] = row[4][0]
                                    if (len(row[3]) != 3):
                                        row[3] = row[3][-3:]
                                    atom_info.append(row)
                        for i in range(len(atom_info)):
                            if (atom_info[i][4][0] == seq):
                                first = i
                                break
                        # print(first)
                        num = 0
                        for i in range(len(atom_info)):
                            if (atom_info[i][4][0] == seq and atom_info[i][0] == 'ATOM'):
                                num += 1
                                last = i
                        # Number of atoms
                        if ((last - first + 1) != num): print(file_name + '-' + seq + '-' + 'error')
                        ff.write(str(last - first + 1))
                        ff.write('\t')
                        ff.write('\n')

                        matrix = np.zeros((last - first + 1, last - first + 1), dtype='float32')

                        mol1 = mol
                        mol2 = mol
                        N1 = mol1.GetNumAtoms()
                        N2 = mol2.GetNumAtoms()
                        xyzs1 = mol1.GetConformer(0).GetPositions()
                        xyzs2 = mol2.GetConformer(0).GetPositions()
                        dismatrix = np.zeros((N1, N2), dtype=np.float16)
                        for i in range(N1):
                            cs = np.tile(xyzs1[i], N2).reshape((N2, 3))
                            dismatrix[i] = np.linalg.norm(xyzs2 - cs, axis=1)

                        xyz = mol.GetConformer(0).GetPositions()

                        # 每个原子
                        for i in range(first, last + 1):
                            degree = 0
                            for j in range(i + 1, last + 1):
                                bond = mol.GetBondBetweenAtoms(i, j)
                                if bond:
                                    matrix[i - first, j - first] = matrix[j - first, i - first] = 1
                                else:
                                    matrix[i - first, j - first] = matrix[j - first, i - first] = 0

                            # print(atom_info[i])
                            if (atom_info[i][4][0] == seq):
                                if n != atom_info[i][5]:
                                    p += 1  # residues num
                                    label = fasta[p]
                                    try:
                                        surface_atom = surface[p]
                                    except:
                                        print()
                                    n = atom_info[i][5]
                                else:
                                    label = fasta[p]
                                    surface_atom = surface[p]

                                # 所在氨基酸编号
                                ff.write(str(p))
                                ff.write('\t')
                                # 标签
                                ff.write(str(label))
                                ff.write('\t')
                                # surface
                                ff.write(str(surface_atom))
                                ff.write('\t')

                                # 原子坐标
                                position = xyz[i].tolist()
                                for m in range(3):
                                    ff.write(str(position[m]))
                                    ff.write('\t')
                                # 原子特征
                                atom = mol.GetAtomWithIdx(i)
                                # print(onek_encoding_unk(atom.GetSymbol(), ELEM_LIST))
                                atom_feature = onek_encoding_unk(atom.GetSymbol(), ELEM_LIST) + onek_encoding_unk(
                                    atom.GetDegree(), [0, 1, 2, 3, 4, 5]) + onek_encoding_unk(atom.GetFormalCharge(),
                                                                                              [-1, -2, 1, 2,
                                                                                               0]) + onek_encoding_unk(
                                    int(atom.GetChiralTag()), [0, 1, 2, 3]) + [1 if atom.GetIsAromatic() else 0]
                                # valence = atom.GetTotalValence()
                                # print(valence)
                                for x in range(39):
                                    if (atom_feature[x] == False):
                                        ff.write('0')
                                        ff.write('\t')
                                    else:
                                        ff.write('1')
                                        ff.write('\t')
                                ff.write('\n')

                        # f1.write(file_name + '\n')
                        point = point + 2
                    else:
                        print(file_name)
                except:
                    print(file_name + '\n')
    print("done extracting feature no pssm")


if __name__ == '__main__':
    main()
