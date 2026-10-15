import sys
import pandas as pd
import os

muts = sys.argv[1]
faspr_path = sys.argv[2]

#XXX: python 4_FASPR_automation_v1.py "./Sample_output/kcat_FASPR_input_v1.csv" <path to FASPR source code>

data = pd.read_csv(muts, sep="\t", header=0)

for i in range(0, len(data)):
        print(data['affin.pdb_id'][i], data['affin.chain'][i], data.mutation[i])
        v1 = data['affin.pdb_id'][i]
        v2 = "./data/PDB_files/"
        
        if(not os.path.exists(v2+v1+".pdb")):
        	v1 = str.lower(v1)
        
        v4 = "./Sample_output/FASPR_input/"
        mutation_list = str.split(data.mutation[i],'/')
        v3=data['affin.chain'][i] + ":" + mutation_list[0]
        for j in range(1, len(mutation_list)):
        	v3 = v3 + "+" + data['affin.chain'][i] + ':' + mutation_list[j]
        print(v1, v2, v3, v4)  #1R37 ../Km_for_mutation/ A:W95L+A:N249Y ./FASPR_IP_Km_mutant/
        os.system('python3 separate_chains_and_seqs_v2.py '+v1+' '+v2+' '+v3+' '+v4)
        os.system(faspr_path+' -i ./Sample_output/FASPR_input/'+v1+'_'+data['affin.chain'][i]+'.pdb -s ./Sample_output/FASPR_input/'+v1+'_'+data['affin.chain'][i]+'_'+v3+'_mut.fasta -o ./Sample_output/FASPR_output/'+v1+'_'+data['affin.chain'][i]+'_'+v3+'_mut.pdb')
        #break

