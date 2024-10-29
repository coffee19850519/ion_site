from fastapi import FastAPI, Form
import subprocess
import os
import shutil
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

# # Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allow frontend origin
    allow_credentials=True,
    allow_methods=["*"],  # Allow all methods (POST, GET, etc.)
    allow_headers=["*"],  # Allow all headers
)

# Test: 
# /api/predict?userId=2024-07-31T16:48:32.585Z5zolqz1qjlk&time=2024-09-10T21:55:20.899Z&pdbId=1A5T_A&ion=ZN&cutoff=0.5
@app.post("/ionpred_api/predict")
async def get_predict_result(userId:str = Form(...), time:str = Form(...), pdbId:str = Form(...), 
                       ion:str = Form(...), cutoff:str = Form(...)):
    print(userId)
    print(time)
    print(pdbId)
    print(cutoff)
    print(ion)

    output_path = f'files/{userId}/{time}/'
    os.makedirs(output_path, exist_ok=True)  

    datalist_path = f'{output_path}data_list.txt'
    pdb_path = f'{output_path}pdb/'
    fasta_path = f'{output_path}fasta/'
    one_pdb_path = f'{output_path}one_pdb/'
    lable_path = f'{output_path}label/'
    pos_file_path = f'{lable_path}pos_file.txt'
    surface_path = f'{output_path}surface/'

    feature_path = f'{output_path}feature/'
    feature_pssm_types_path = f'{feature_path}feature_pssm_types.txt'
    feature_pssm_path = f'{feature_path}feature_pssm.txt'
    feature_no_pssm_path = f'{feature_path}feature_no_pssm.txt'

    pssm_path = f'{output_path}pssm/'

    result_path = f'{output_path}result/'
    api_output_path = f'{output_path}api_output.txt'

    scripts_with_params = [
        ('get_data_list.py', ['-pdb', pdbId, '-opath', datalist_path]),
        ('get_pdb.py', ['-input', datalist_path, '-output', pdb_path]),
        ('divide_pdb.py', ['-input', datalist_path, '-pdb-path', pdb_path,'-one-pdb-path', one_pdb_path]),
        ('get_fasta.py', ['-input', datalist_path, '-one-pdb-path', one_pdb_path, '-fasta-dir', fasta_path]),
        ('new_positive.py', ['-ion', ion, '-input', datalist_path, '-label-path', lable_path, '-pos-file-path', pos_file_path, '-pdb-path', pdb_path]),
        ('get_label.py', ['-input', datalist_path, '-label-path', lable_path, '-one-pdb-path', one_pdb_path, '-pos-file-path', pos_file_path, '-ion', ion]),
        ('get_surface.py', ['-input', datalist_path, '-surface-file', surface_path,'-one-pdb-path', one_pdb_path, '-ion', ion]),
        ('get_pssm.py', ['-input', datalist_path, '-blast', './tool/ncbi_blast/blast-BLAST_VERSION+/db/', '-in-pssm', fasta_path, '-out-pssm', pssm_path]),
        ('get_feature_no_pssm.py', ['-ipath', output_path, '-label', lable_path, '-surface', surface_path, '-feature', feature_path, '-one-pdb-path', one_pdb_path,]),
        ('add_pssm.py', ['-ipath', output_path, '-label', lable_path, '-pssm-path', pssm_path, '-feature-pssm', feature_pssm_path, '-feature-no-pssm', feature_no_pssm_path]),
        ('add_types.py', ['-fasta', fasta_path, '-cat-feature', feature_pssm_path, '-feature-pssm-types-path', feature_pssm_types_path, '-feature', feature_path, '-ion', ion]),
        ('test.py', ['-ipath', feature_pssm_types_path, '-ion', ion, '-result-path', result_path]),
        ('get_result_output.py', ['-datalist', datalist_path, '-ion', ion, '-Cutoff', '0.5', '-result-path', result_path, '-onepdb-path', one_pdb_path,'-api-output-path', api_output_path,]),
    ]

    try:
        for script, params in scripts_with_params:
            run_script(script, params)
      
        with open(api_output_path, 'r') as file:
            content = file.read()
        
        shutil.rmtree(output_path)
        # Return the content after successful cleanup
        return content
            
    except FileNotFoundError:
        return "File not found."
    except Exception as e:
        return f"An error occurred: {e}"
     

# 定义要执行的脚本列表和顺序，以及每个脚本的参数

def run_script(script_name, params):
    try:
        python_executable = 'python3'
        # 构建执行命令，包括参数
        command = [python_executable, script_name] + params
        # 执行脚本
        subprocess.run(command, check=True)
        # print(f"{script_name} executed successfully with parameters {params}.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing {script_name}: {e}")