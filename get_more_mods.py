import os
import glob
import argparse
import shutil

'''
Search out one_pdb files with multiple mods
'''
def main():
    parser = argparse.ArgumentParser(description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-ipath', dest='ipath', type=str, help='Specify the directory to search', required=True)
    parser.add_argument('-opath', dest='opath', type=str, help='Specify a new file for the one_pdb name with multiple mods', required=True)
    parser.add_argument('-npath', dest='npath', type=str,help='Specify the path where the files of multiple mods are newly stored', required=True)
    args = parser.parse_args()


    ipath = args.ipath
    opath = args.opath
    npath = args.npath
    # os.makedirs(ipath, exist_ok=True)
    os.makedirs(npath, exist_ok=True)

    fw = open(opath, 'w')
    ipath = ipath
    for file_abs in glob.glob(ipath):
        # print(file_abs)
        # list.append(i)
        # print(file_abs)
        max = 0

        with open(file_abs,'r')as f:
            for line in f:
                row = line.strip().split()
                if len(row) == 2 and row[0] == 'MODEL':
                    if max < int(row[1]):
                        max = int(row[1])
        if max != 0:
            file_path, file_name = os.path.split(file_abs)
            # os.rename(file_name,file_name[:-3]+ '-mods' + '.one_pdb')
            if not os.path.exists(npath):
                os.makedirs(npath)
            shutil.move(file_abs, npath + file_name)
            print("move %s -> %s" % (file_abs, npath + file_name))

            # print(file_abs[-10:-4] + ':' + str(max) + '\n')
            fw.write(file_abs[-10:-4] + ':' + str(max) + '\n')
    print("done generating more models!")
if __name__ == "__main__":
    main() 
    

