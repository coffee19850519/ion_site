import requests
import os
import argparse
from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser

# Download the PDB files from the RCSB database website

# def remove_pdb_null(filePath):
#   allfiles = os.listdir(filePath)
#   fw = open('./data/null_pdb.txt','w')
#   for file in allfiles:
#       itemPath = os.path.join(filePath, file)
#       if not os.path.isdir(itemPath):
#           # 获取文件的大小
#           fileSize = os.path.getsize(itemPath)
#           if fileSize <= 260:
#               # print(f'该文件的大小为{fileSize}字节，路径为{itemPath}')
#               fw.write(str(itemPath.strip().split('.')[1][-4:]) + '\n')
#       else:
#           remove_pdb_null(itemPath)

def main():
  parser=argparse.ArgumentParser(description='supply the input and outpaths to download the one_pdb files for a specific ion i.e. ZN, CA, CO3')
  parser.add_argument('-input', dest='datalist_path', type=str, help='Specify the location of one_pdb list i.e. data_list.txt.txt file for the ion of interest', required=True)
  parser.add_argument('-output', dest='output_path', type=str, help='Specify the location for the one_pdb files to be stored', required=True)

  args = parser.parse_args()

  datalist_path = args.datalist_path
  output_path = args.output_path
  os.makedirs(output_path, exist_ok=True)

  with open(datalist_path, 'r') as fp:
    lines = fp.read().split('\n')
    pre_id = ''
    for line in lines:
      id = line[:4]
      if id != pre_id:
        pre_id = id
        try:
          r = requests.get('https://files.rcsb.org/download/' + str(id) + '.pdb', stream=True)
          print(format(id), "already downloaded", )
          with open(output_path + str(id) + '.pdb', 'wb') as f:
            f.write(r.content)
        except:
          print('could not download pdb for '+id)
          # fw.write(str(id) + '\n')
  print("Done downloading PDB files!")

  # remove_pdb_null(output_path[:-1])

if __name__ == '__main__':
  main()