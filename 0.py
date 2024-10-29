import sys
# print("Python executable:", sys.executable)
'''
参数：PDB文件路径、输出文件路径、离子类型、cutoff
'''
import argparse
import subprocess

parser = argparse.ArgumentParser(
    description='')
parser.add_argument('-pdb_path', dest='pdb_path', type=str, help='PDB file path', default='./data/pdb/1A5T_A.pdb')

parser.add_argument('-output_path', dest='output_path', type=str, help='Output file path',default='./data/output/Prediction_results.txt')

parser.add_argument('-ionic_type', dest='ionic_type', type=str, help='Ionic type', default='ZN')

parser.add_argument('-cutoff', dest='cutoff', type=str, help='cutoff',default='0.5')

args = parser.parse_args()

pdb_chain_path = args.pdb_path
output_path = args.output_path
ion = args.ionic_type
cutoff = args.cutoff

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

# 定义要执行的脚本列表和顺序，以及每个脚本的参数
scripts_with_params = [
    ('get_data_list.py', ['-pdb', pdb_chain_path]),
    ('get_pdb.py', ['-input', 'data_list.txt', '-output', 'data/'+ion+'/pdb/']),
    ('get_fasta.py', ['-input', 'data_list.txt', '-fasta-dir', 'data/'+ion+'/fasta/', '-opath', 'data/'+ion+'/one_pdb/']),
    ('new_positive.py', ['-ion', ion]),
    ('get_label.py', ['-input', 'data_list.txt', '-label-path', 'data/'+ion+'/label/', '-ion', ion]),
    ('get_surface.py', ['-input', 'data_list.txt', '-surface_file', 'data/'+ion+'/surface/', '-ion', ion]),
    ('get_pssm.py', ['-input', 'data_list.txt', '-blast', './tool/ncbi_blast/blast-BLAST_VERSION+/db/', '-in-pssm', 'data/'+ion+'/fasta/', '-out-pssm', 'data/'+ion+'/pssm/']),
    ('get_feature_no_pssm.py', ['-ipath', 'data/'+ion+'/', '-label', 'data/'+ion+'/label/', '-surface', 'data/'+ion+'/surface/', '-feature', 'data/'+ion+'/feature/']),
    ('add_pssm.py', ['-ipath', 'data/'+ion+'/', '-label', 'data/'+ion+'/label/', '-feature', 'data/'+ion+'/feature/']),
    ('add_types.py', ['-fasta', 'data/'+ion+'/fasta/', '-cat-feature', 'data/'+ion+'/feature/feature_pssm.txt', '-feature', 'data/'+ion+'/feature/', '-ion', ion]),
    ('test.py', ['-ipath', 'data/'+ion+'/feature/feature_pssm_types.txt', '-ion', ion]),
    ('get_result_output.py', ['-datalist', 'data_list.txt', '-ion', ion, '-Cutoff', '0.5']),
]

def run_script(script_name, params):
    try:
        python_executable = 'python3'
        # 构建执行命令，包括参数
        command = [python_executable, script_name] + params
        # 执行脚本
        subprocess.run(command, check=True)
        print(f"{script_name} executed successfully with parameters {params}.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing {script_name}: {e}")

def main():
    for script, params in scripts_with_params:
        run_script(script, params)

if __name__ == '__main__':
    main()