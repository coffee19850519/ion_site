import argparse

parser = argparse.ArgumentParser(
    description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')

parser.add_argument('-datalist', dest='datalist', type=str, help='Specify the path to  the specific ion of interest',
                    required=True)
parser.add_argument('-ion', dest='ion', type=str, help='Specify the directory for the ion under consideration',
                    required=True)
parser.add_argument('-Cutoff', dest='Cutoff', type=str, help='Specify the path to the features',
                    required=True)
parser.add_argument('-result-path', dest='result_path', type=str, help='Specify the path to  the specific ion of interest',
                    required=True)
parser.add_argument('-onepdb-path', dest='onepdb_path', type=str, help='Specify the path to  the specific ion of interest',
                    required=True)
parser.add_argument('-api-output-path', dest='api_output_path', type=str, help='Specify the path to  the specific ion of interest',
                    required=True)

args = parser.parse_args()

data_list = args.datalist
ion = args.ion
Cutoff = args.Cutoff
result_path = args.result_path
onepdb_path = args.onepdb_path
api_output_path = args.api_output_path

with open(data_list,'r')as f:
    lines=f.readlines()
data=lines[0].split()[0]
# data='1A5T_A'
# data='>'+ data.replace('_','-')
chain=data[-1]

# ion = 'ZN'
#
Cutoff=float(Cutoff)

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

residues_name= {'ZN': 'Zinc',
                'CU': 'Copper',
                'FE2': 'Ferrous ',
                'CA': 'Calcium',
                'MG': 'Magnesium',
                'MN': 'Manganese',
                'NA': 'Sodium',
                'K': 'Potassium'}

candidate_num=0
fw=open(api_output_path, 'w')
fw.write('ID'+'\t'+'Position'+'\t'+'Residue'+'\t'+'scores'+'\t'+f'Cutoff={Cutoff}'+'\t'+'\n')

with open(f'{result_path}result_0.txt', 'r')as fs:
    lines=fs.readlines()

sequence = []
data_list = []
with open(f'{onepdb_path}{data}.pdb', 'r') as f:
    m = -10000
    for line1 in f:
        row = line1.strip().split()
        if (row[0] == 'ATOM' and row[4][0] == chain):
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

            if str(m) != row[5]:  # 氨基酸改变
                if row[3] == 'GLY':
                    res='G'
                elif row[3] == 'ALA':
                    res='A'
                elif row[3] == 'VAL':
                    res='V'
                elif row[3] == 'LEU':
                    res='L'
                elif row[3] == 'ILE':
                    res='I'
                elif row[3] == 'PRO':
                    res='P'
                elif row[3] == 'PHE':
                    res='F'
                elif row[3] == 'TYR':
                    res='Y'
                elif row[3] == 'TRP':
                    res='W'
                elif row[3] == 'SER':
                    res='S'
                elif row[3] == 'THR':
                    res='T'
                elif row[3] == 'CYS':
                    res='C'
                elif row[3] == 'MET':
                    res='M'
                elif row[3] == 'ASN':
                    res='N'
                elif row[3] == 'GLN':
                    res='Q'
                elif row[3] == 'ASP':
                    res='D'
                elif row[3] == 'GLU':
                    res = 'E'
                elif row[3] == 'LYS':
                    res='K'
                elif row[3] == 'ARG':
                    res='R'
                elif row[3] == 'HIS':
                    res='H'
                else:
                    res='X'
                m = str(row[5])
                pos = row[5]
                sequence.append(str(res))
                score = 0
                if res in candidate_residues[ion]:
                    # pos = row[5]
                    score = lines[candidate_num].split()[-1][:5]
                    candidate_num+=1

                entry = {
                    "id": data,
                    "position": str(pos),
                    "residue": str(res),
                    "score": str(score)
                }
                data_list.append(entry)
                
fw.write('>'+ data + '\t' +''.join(sequence) +'\n')
name = residues_name.get(ion)
for entry in data_list:
    line = f"{entry['id']}\t{entry['position']}\t{entry['residue']}\t{name}:{entry['score']}\n"
    fw.write(line)
