import argparse
import os

parser = argparse.ArgumentParser(
        description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
parser.add_argument('-pdb_chain', dest='pdb_chain', type=str, help='',
                    required=True)
parser.add_argument('-opath', dest='opath', type=str, help='Specify the path to save the datalist.txt',
                        required=True)

args = parser.parse_args()

pdb_chain_path = args.pdb_chain
output_path = args.opath
pdb_chain = pdb_chain_path.split('/')[-1].split('.')[0]

# os.makedirs('data/ZN/', exist_ok=True)

fw = open(output_path,'w')
fw.write(pdb_chain)
fw.close()