from biopandas.pdb import PandasPdb
import os
import method as md
cwd = os.getcwd()
ppdb = PandasPdb()

#main takes in input a folder, which contains a subfolder for each protein that i want to compute. 
#The subfolder is named as the name of the protein and contains the file PDB AlphaFold2 and a folder named "wt" which contains the experimental file PDB.
def main(proteins_folder, dist=5):
        
    #proteins_folder contiene una cartella per ogni proteina considerata
    proteins = os.listdir(proteins_folder)

    #scorro tutte le proteine
    for p_name in proteins:
        print(p_name)
        path_af = proteins_folder + '\\' + p_name + '\\' + p_name + '_af.pdb'
        mut = proteins_folder+ '\\' + p_name + '\\mut'
        path_exp = mut + '\\' + os.listdir(mut)[0]
        path_alignment = proteins_folder + '\\' + p_name + '\\' + p_name + '_alignment.pdb'
        path_out = proteins_folder + '\\' + p_name + '\\' + p_name + '_out.txt'
        
        md.apply(p_name, path_af, path_exp, path_alignment, path_out, dist)
        
    return   
