#Questo file contiene le funzioni per il secondo e per il terzo step del metodo, ovvero i calcoli su tutta la proteina e sul sito attivo.

from biopandas.pdb import PandasPdb
import pandas as pd
import os
import numpy as np
cwd = os.getcwd()
ppdb = PandasPdb()

#rmsd calcola l'rmsd tra due dataframe
def rmsd(df_af, df_exp):
    df_af = df_af.rename(columns={"x_coord_af": "x_coord", 'y_coord_af' : 'y_coord', 'z_coord_af':'z_coord'})
    df_exp = df_exp.rename(columns={"x_coord_exp": "x_coord", 'y_coord_exp' : 'y_coord', 'z_coord_exp':'z_coord'})
    if (df_af.shape[0]==0 & df_exp.shape[0]==0):
        return 0
    else:
        return ppdb.rmsd(df_af, df_exp)
    
def get_sideChains(df):
    #lista di atomi della backbone
    bb_lis = ['CA','HA','N','C','O','HN','H']   
    
    #le righe di ser contangono true solo se l'atomo appartiene a una catena laterale
    ser1 = (~df.iloc[:,3].isin(bb_lis))
    #side_chains contiene solo gli atomi del dataframe appartenenti a catene laterali, ma anche gli idrogeni
    side_chainsT = df.loc[ser1]
    
    #le righe di contengono true solo se l'atomo non è un idrogeno
    ser2= ~side_chainsT.iloc[:,3].str.startswith("H")
    #side_chains ontiene solo gli atomi del dataframe appartenenti a catene laterali, senza gli idrogeni
    side_chains = side_chainsT.loc[ser2]
    return side_chains

def get_backbone(df):
    #lista di atomi della backbone
    bb_lis = ['CA','HA','N','C','O','HN','H']   
    
    #le righe di ser contangono true solo se l'atomo appartiene a una backbone
    ser1 = (df.iloc[:,3].isin(bb_lis))
    #side_chains contiene solo gli atomi del dataframe appartenenti alla backbone, ma anche gli idrogeni
    backbone_T = df.loc[ser1]
    
    #le righe di contengono true solo se l'atomo non è un idrogeno
    ser2= ~backbone_T.iloc[:,3].str.startswith("H")
    #bakbone ontiene solo gli atomi del dataframe appartenenti a catene laterali, senza gli idrogeni
    backbone = backbone_T.loc[ser2]
    return backbone

def calculate_differences(join):

    all_residues_af = list(join['residue_number_af'].drop_duplicates())
    
    dict_amino = {}
        
    for res in all_residues_af:      
        #join considerando solo il residuo ennesimo
        join_res = join.loc[join.residue_number_af==res]
        #standard deviation dell'intero residuo
        join_res_std = join_res['distance'].std() 
        
        #dataframe 
        df_af_res = join_res.loc[: , 'record_name_af':'line_idx_af']
        df_exp_res = join_res.loc[: , 'record_name_exp':'line_idx_exp']
        
        #sidechains
        df_af_res_sc = get_sideChains(df_af_res)
        df_exp_res_sc = get_sideChains(df_exp_res)
        
        join_res_SC = get_sideChains(join_res) #sidechains del residuo
           
        len_SC = len(join_res_SC)
        join_res_SC_std = join_res_SC['distance'].std() #standard deviation delle sidechains del residuo

        dict_amino[str(res)] = [df_af_res.iloc[0]['residue_name_af'], df_exp_res.iloc[0]['residue_name_exp'], rmsd(df_af_res, df_exp_res), rmsd(df_af_res_sc, df_exp_res_sc), join_res_std, join_res_SC_std, len_SC] 
    
    df_diff = pd.DataFrame.from_dict(dict_amino, orient="index").rename(columns={0:'RESIDUE_NAME_AF', 1:'RESIDUE_NAME_EXP', 2:'TOTAL_RMSD', 3:'SIDECHAIN_RMSD', 4:'TOTAL_STD', 5:'SIDECHAIN_STD', 6:'SIDECHAIN_ATOMS'})
    df_diff.index.name = "residue_number_af"
    
    return df_diff

def zone(path_alignment, path_exp, dist):
    ppdb.read_pdb(path_exp)
    df_ligand = ppdb.df['HETATM']
    arr_res = []
    
    #elimino le molecole d'acqua
    df_ligand=df_ligand.loc[df_ligand.residue_name!='HOH']    
    
    #scorro gli atomi del ligando
    for l in df_ligand.index:
        raw_l = df_ligand.iloc[l] 
        
        #coordinate dell'atomo
        l_coordinates = (raw_l.loc['x_coord'],raw_l.loc['y_coord'],raw_l.loc['z_coord'] )
        
        #calcolo le distanze dal modello di alfafold2 dall'atomo
        ppdb.read_pdb(path_alignment)
        distances = ppdb.distance(xyz = l_coordinates, records=('ATOM'))
        
        #filtro gli atomi con distanza < dist e aggiungo i residui che hanno almeno un atomo all'interno della zona alla lista
        all_within_dist = ppdb.df['ATOM'][distances < dist]
        arr_res += all_within_dist['residue_number'].drop_duplicates().tolist()
        
    arr_res = np.unique(arr_res) 
    return arr_res
