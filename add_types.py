
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-fasta', dest='fasta', type=str, help='Specify the path to  the specific ion of interest', required=True)
    parser.add_argument('-cat-feature', dest='cfeature', type=str, help='Specify the directory for the ion under consideration', required=True)
    parser.add_argument('-feature-pssm-types-path', dest='feature_pssm_types_path', type=str, help='Specify the directory for the ion under consideration', required=True)
    parser.add_argument('-feature', dest='feature', type=str, help='Specify the path to the features', required=True)
    parser.add_argument('-ion', dest='ion', type=str, help='Specify the path to the features', required=True)
    
    args = parser.parse_args()

    fasta = args.fasta
    cfeature = args.cfeature
    feature = args.feature
    feature_pssm_types_path = args.feature_pssm_types_path
    ion = args.ion
    # res_dict={'ZN':'C,D,H,E'}
    res_dict={'ZN': 'C,D,E,H',
              'CU': 'C,H',
              'FE2': 'D,E,H',
              'CA': 'D,E,G,N',
              'MG': 'D,E,N',
              'MN': 'D,E,H',
              'NA': 'D,G,N,S',
              'K': 'D,E,G,N,S'}
    
    residues = list(res_dict[ion].replace(",",''))

    num = 0
    # wn = open('./hh_why_colabfold.txt','w')
    with open(fasta + 'fasta.txt', 'r') as la:
        with open(feature_pssm_types_path, 'w') as nfe:
            with open(cfeature, 'r') as fe:
                num_acid = fe.readline().strip()
                nfe.write(num_acid + '\n')
                for i in range(int(num_acid)):
                    name_pssm = la.readline().strip().split()[0]
                    la.readline()
                    name_lines = name_pssm[1:].split('-')
                    # name = name_lines[0] + '-' + name_lines[1]
                    name = name_lines[0]
                    chain = name_lines[1]
                    num_fa = 0

                    with open(fasta + str(name) + '-' + str(chain) + '.fa', 'r') as ps:
                        ps.readline()
                        # ps.readline()
                        # ps.readline()
                        num_atom = fe.readline().strip()
                        nfe.write(num_atom + '\n')
                        line = ps.readline().strip().split()

                        # if name == '4KO1':
                        #     print()
                        for i_a in range(int(num_atom.split()[0])):
                            # print(line[0])
                            # try:
                            #      # seq_ac_pssm = int(line[0]) - 1
                            #      num_fa = num_fa - 1
                            # except:
                            #     print(name)
                            feature_a = fe.readline().strip()
                            # if name == '4KO1':
                            #     wn.write(feature_a + '\t')
                            #     wn.write('\n')
                            seq_ac = feature_a.split()[0]
                            if num_fa != int(seq_ac):
                                # line = ps.readline().strip().split()
                                num_fa += 1
                            # if name == '4KO1' and num_fa == 486:
                            #     print()
                            # print(num_fa)
                            try:# 283
                                fa = line[0][num_fa]
                            except:
                                print(line[0])
                            # print(fa)
                            if fa in residues:
                                fa = '1'
                            else:
                                fa = '0'
                            nfe.write(fa + '\t')
                            nfe.write(feature_a + '\t')
                            # for num_p in range(pssm.__len__()):
                            #     nfe.write(pssm[num_p] + '\t')
                            nfe.write('\n')

                            # print(str(name) + '-' + str(chain))
                            # print(num_fa)



    print("Done adding types!")   
        
if __name__ == "__main__":
    main() 
