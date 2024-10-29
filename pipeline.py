import os
# os.system('python get_pdb.py -input data/CU/data_list.txt -output data/CU/one_pdb/')
'''
 通过网页用户选择，得到参数 ion、PDB文件、pdb_chain格式的输入
'''
os.system('python get_data_list.py -input pdb_chain')
os.system('python pdb_update.py -input data/data_list.txt -opath data/pdb/ -npath data/pdb_update/')
os.system('python divide_pdb.py -input data/data_list.txt -opath data/pdb_update/ -npath data/one_pdb/')
os.system('python get_more_mods.py -ipath data/one_pdb/* -opath data/more_mods.txt -npath data/more_mods_pdb/')
os.system('python remove_mods.py -input data/more_mods.txt -mpath data/more_mods_pdb/ -output data/one_pdb/')
os.system('python get_fasta.py -input data/data_list.txt -fasta-dir data/fasta/ -model-one-path data/one_pdb/')
os.system('python get_label.py -input data/data_list.txt -label-path data/label/label.txt')
os.system('python get_surface.py -input data/data_list.txt -surface-dir data/surface/surface.txt')
os.system(r'python get_pssm.py -input data/data_list.txt -blast D:\LenovoSoftstore\Install\ncbi_blast\blast-BLAST_VERSION+\db\ -in-pssm data/fasta// -out-pssm data/pssm//')
os.system('python get_feature_no_pssm.py -ipath data/ -label data/label/ -surface data/surface/ -feature data/feature/')
os.system('python add_pssm.py -ipath data/ -label data/label/ -feature data/feature/ ')
os.system('python add_types.py -fasta data/fasta/ -cat-feature data/feature/feature_pssm.txt  -feature data/feature/ -ion ZN')
os.system('python test.py -input data/feature/feature_pssm_types.txt -ion ZN')
