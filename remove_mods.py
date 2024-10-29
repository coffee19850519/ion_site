import argparse
import os
def remove_more_mods(input,mpath, output):
    with open(input,'r')as f:
        for line in f:
            row = line.strip().split(':')
            pdb = row[0]
            fw = open(output+ str(pdb) + '.pdb','w')
            with open(mpath+ str(pdb) + '.pdb','r')as fp:
                li = fp.readlines()
            for count in range(len(li)):
                if li[count].strip().split()[0] == 'MODEL' and li[count].strip().split()[1] == '2':
                    break
                fw.write(li[count].strip().split('\n')[0] + '\n')
                
                        
def main():
    parser = argparse.ArgumentParser(description='supply the old and new directory to update the downloaded pdbs for a specific ion i.e. ZN, CA, CO3')
    parser.add_argument('-input', dest='input', type=str, help='more_model_total file', required=True)
    parser.add_argument('-mpath', dest='mpath', type=str, help='Specify the location of one_pdb list of more mods i.e. more_mods.txt file for the ion of interest', required=True)
    parser.add_argument('-output', dest='output', type=str, help='Specify the location for the output', required=True)
    
    args = parser.parse_args()
    
    
    input = args.input
    mpath = args.mpath
    output = args.output
    
    os.makedirs(output, exist_ok=True)
    
    remove_more_mods(input, mpath, output)
    print("Done removing model!")
    

if __name__ == '__main__':
  main()
  
       
