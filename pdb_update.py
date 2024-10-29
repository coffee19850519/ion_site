# Remove the 'Name of non-standard residue' line from the downloaded one_pdb file, starting each line with 'HETATM'
import requests
import os
import argparse
from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser

def remove_HETATM(opath, npath):
    with open(opath, 'r') as op:
        write_all = ''
        for line in op:
            if line[:6] != 'HETATM':
                write_all += line  
    with open(npath, 'w') as np:
        np.write(write_all)       
    
    
def main():
    parser = argparse.ArgumentParser(description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-input', dest='input', type=str, help='Specify the location of file that contain the one_pdb chains i.e. data_list.txt.txt for the specific ion of interest', required=True)
    parser.add_argument('-opath', dest='opath', type=str, help='Specify the location of one_pdb list i.e. data_list.txt.txt file for the ion of interest', required=True)
    parser.add_argument('-npath', dest='npath', type=str, help='Specify the location for the one_pdb files to be stored', required=True)

    args = parser.parse_args()

    input = args.input
    opath = args.opath
    npath = args.npath
    os.makedirs(npath, exist_ok=True)
    
    with open(input, 'r') as f:
        lines = f.read().split()
        for line in lines:
            name = line[:4].upper()
            # print(name)
            remove_HETATM(opath + str(name) + '.pdb', npath + str(name) + '.pdb')
    print("completed pdb_update!")

if __name__ == "__main__":
    main() 


