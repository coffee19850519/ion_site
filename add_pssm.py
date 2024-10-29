import os
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-ipath', dest='ipath', type=str, help='Specify the path to  the specific ion of interest',
                        required=True)
    parser.add_argument('-label', dest='label', type=str, help='Specify the directory for the ion under consideration',
                        required=True)
    parser.add_argument('-pssm-path', dest='pssm_path', type=str, help='Specify the path to the features', required=True)
    parser.add_argument('-feature-pssm', dest='feature_pssm', type=str, help='Specify the path to the features', required=True)
    parser.add_argument('-feature-no-pssm', dest='feature_no_pssm', type=str, help='Specify the path to the features', required=True)

    args = parser.parse_args()

    ipath = args.ipath
    label = args.label
    feature_pssm = args.feature_pssm
    feature_no_pssm = args.feature_no_pssm
    pssm_path = args.pssm_path
    
    os.makedirs(pssm_path, exist_ok=True)

    num = 0
    with open(label + 'label.txt', 'r') as la:
        with open(feature_pssm, 'w') as nfe:
            with open(feature_no_pssm, 'r') as fe:
                num_acid = fe.readline().strip()
                nfe.write(num_acid + '\n')
                for i in range(int(num_acid)):
                    name_pssm = la.readline().strip().split()[0]
                    la.readline()
                    name_lines = name_pssm[1:].split('_')
                    # name = name_lines[0] + '-' + name_lines[1]
                    name = name_lines[0]

                    chain = name_lines[1]

                    with open(pssm_path + str(name) + '-' + str(chain) + '.pssm', 'r') as ps:
                        ps.readline()
                        ps.readline()
                        ps.readline()
                        num_atom = fe.readline().strip()
                        nfe.write(num_atom + '\n')
                        line = ps.readline().strip().split()
                        for i_a in range(int(num_atom.split()[0])):
                            # print(line[0])

                            seq_ac_pssm = int(line[0]) - 1

                            # print(name)

                            feature_a = fe.readline().strip()
                            nfe.write(feature_a + '\t')
                            seq_ac = feature_a.split()[0]
                            if seq_ac_pssm != int(seq_ac):
                                line = ps.readline().strip().split()
                                seq_ac_pssm += 1
                            pssm = line[2:22]
                            for num_p in range(pssm.__len__()):
                                nfe.write(pssm[num_p] + '\t')
                            nfe.write('\n')

    print("Done adding PSSM!")


if __name__ == "__main__":
    main() 