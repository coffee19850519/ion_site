import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-input', dest='input', type=str, help='Specify the location of file that contain the one_pdb chains i.e. data_list.txt.txt for the specific ion of interest', required=True)
    parser.add_argument('-label-path', dest='lpath', type=str, help='Specify the directory for the ion under consideration', required=True)
    parser.add_argument('-one-pdb-path', dest='one_pdb_path', type=str, help='Specify the directory for the ion under consideration', required=True)
    parser.add_argument('-pos-file-path', dest='pos_file_path', type=str, help='Specify the path to the features', required=True)
    parser.add_argument('-ion', dest='ion', type=str, help='Specify the file containing positive chains for the ion under consideration', required=True)
    
    args = parser.parse_args()

    input = args.input
    lpath = args.lpath
    one_pdb_path = args.one_pdb_path
    pos_file_path = args.pos_file_path
    ion = args.ion

    os.makedirs(lpath, exist_ok=True)

    pos_file = pos_file_path
    
    index = input.rfind('/')
    one_pdb = input[:index]
     
    fw = open(f'{lpath}label.txt', 'w')
    with open(input, 'r')as f:
        for li in f:
            name = li.strip().split()[0][:4]
            chain = li.strip().split()[0][5]
            # residue = li.strip().split()[1]
            # index = li.strip().split()[3]
            label = []
            with open(f'{one_pdb_path}{str(name)}_{str(chain)}.pdb','r')as fp:
                m = -10000
                for line in fp:
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
                            row[8] = row[7]
                            row[7] = row[6]
                            row[6] = row[5]
                            row[5] = row[4]
                            row[4] = row[3]
                            row[3] = row[2][4:]
                            row[2] = row[2][:4]
                        if (len(row[4]) != 1):
                            row.append('0')
                            # row[11] = row[10]
                            row[7] = row[6]
                            row[6] = row[5]
                            row[5] = row[4][1:]
                            row[4] = row[4][0]
                        if (len(row[3]) != 3):
                            row[3] = row[3][-3:]
                        # fw.write(str(row)+'\n')


                        if str(m) != row[5]:  # 氨基酸改变
                            flag = 0
                            m = str(row[5])
                            if row[3] == 'GLY':row[3] = 'G'
                            elif row[3] == 'ALA':row[3] = 'A'
                            elif row[3] == 'VAL':row[3] = 'V'
                            elif row[3] == 'LEU':row[3] = 'L'
                            elif row[3] == 'ILE':row[3] = 'I'
                            elif row[3] == 'PRO':row[3] = 'P'
                            elif row[3] == 'PHE':row[3] = 'F'
                            elif row[3] == 'TYR':row[3] = 'Y'
                            elif row[3] == 'TRP':row[3] = 'W'
                            elif row[3] == 'SER':row[3] = 'S'
                            elif row[3] == 'THR':row[3] = 'T'
                            elif row[3] == 'CYS':row[3] = 'C'
                            elif row[3] == 'MET':row[3] = 'M'
                            elif row[3] == 'ASN':row[3] = 'N'
                            elif row[3] == 'GLN':row[3] = 'Q'
                            elif row[3] == 'ASP':row[3] = 'D'
                            elif row[3] == 'GLU':row[3] = 'E'
                            elif row[3] == 'LYS':row[3] = 'K'
                            elif row[3] == 'ARG':row[3] = 'R'
                            elif row[3] == 'HIS':row[3] = 'H'

                            # fw.write(str(row)+'\n')

                            with open(pos_file,'r')as ff:
                                for lii in ff:
                                    name1 = lii.strip().split()[0][:4]
                                    chain1 = lii.strip().split()[0][5]
                                    residue = lii.strip().split()[1]
                                    index = lii.strip().split()[3]
                                    if name == name1 and chain == chain1:
                                        if row[5] == index and row[3] == residue:
                                            # label.append('1')
                                            flag = 1
                                            break



                            if flag == 1:
                                label.append('1')
                            else:
                                label.append('0')
                            flag=2
            label = ''.join(str(i) for i in label)
            fw.write('>' + str(name) + '_' + str(chain) + '\n')
            fw.write(str(label) + '\n')
            print(str(name) + 'over' + '\n')
    print("Done extracting label")
if __name__ == '__main__':

    main()

