from biopandas.pdb import PandasPdb
import pandas as pd
import os
import subprocess as sp
import numpy as np
cwd = os.getcwd()
ppdb = PandasPdb()

#eliminateH(df) prende un dataframe che rappresenta una proteina e restituisce il dataframe senza gli atomi di idrgeno
def eliminateH(df):
    #le righe di temp contengono true solo se l'atomo è un idrogeno
    temp = ~df['atom_name'].str.startswith("H")
    ris = df.loc[temp]
    return ris.reset_index(drop=True)

#rmsd calcola l'rmsd tra due dataframe
def rmsd(df_af, df_exp):
    df_af = df_af.rename(columns={"x_coord_af": "x_coord", 'y_coord_af' : 'y_coord', 'z_coord_af':'z_coord'})
    df_exp = df_exp.rename(columns={"x_coord_exp": "x_coord", 'y_coord_exp' : 'y_coord', 'z_coord_exp':'z_coord'})
    if (df_af.shape[0]==0 & df_exp.shape[0]==0):
        return 0
    else:
        return ppdb.rmsd(df_af, df_exp)
    
#distance pende in input due righe che rappresentano atomi e ritorna a loro distanza
def distance(atom1, atom2):
    p1 = np.array([atom1.iloc[0], atom1.iloc[1], atom1.iloc[2]])
    p2 = np.array([atom2.iloc[0], atom2.iloc[1], atom2.iloc[2]])
    #print(p1,p2)
    squared_dist = np.sum((p1-p2)**2, axis=0)
    return np.sqrt(squared_dist)

#get_join(df_exp, df_af) restituisce la join dei due dataframe e aggiunge la colonna 'distance'
def get_join(df_af, df_exp):
    join = df_af.join(df_exp, lsuffix='_af', rsuffix='_exp')
    dist = []
    for i in range(len(join)):
        atom1 = join.iloc[i].loc['x_coord_af':'z_coord_af']        
        atom2 = join.iloc[i].loc['x_coord_exp':'z_coord_exp'] 
        dist.append(distance(atom1, atom2))
    join['distance']= dist
    return join

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

#writeScript compone lo script e lo esegue su ChimeraX. Prende in input il file pdb, il file di AF, il file di output e il file dove scrivere lo script.
def writeScript(path_exp, path_af, path_alignment, path_script):
    
    file_out = open(path_alignment, "w+")
    filecontentX = "open '" + path_exp + "'\n"
    filecontentX += "select /B\ndelete atoms sel\ndelete bonds sel\n"
    filecontentX += "open '" + path_af + "'\n"
    filecontentX += "matchmaker #2 to #1\n"
    filecontentX += "save '" + path_alignment + "' models #2 relModel #1\n"
    filecontentX += "exit"
    file_out.close()
    
    fileX = open(path_script, "w+")
    fileX.write(filecontentX)
    fileX.close()
    
    sp.call(["C:\\Program Files\\ChimeraX\\bin\\ChimeraX.exe", path_script])

#main   
def main(proteins_folder, dist=5):
    
    #proteins_folder contiene una cartella per ogni proteina considerata
    proteins = os.listdir(proteins_folder)
    
    dic = {'3f8y':'Bad','1v2k':'Good','1ql7':'Good','1oh0':'Good','1k3l':'Bad','1nje':'Bad', '2ygf':'Bad', '1syh':'Bad', '1ydk':'*Good','1nja':'Bad','3kgu':'Good','3cwk': 'Good','1a4h':'Bad','2fr3' : '*Bad','3fs6' :'Bad','1njc' :'Bad','2yge': '*Bad','3neo' :'Good','3m4h' :'Bad','3nex':'Good','1syi':'Good','2pdy':'Bad','2fzd':'Bad','2g78':'Good','3kgt':'Good','2pdg':'*Bad','2pdl':'Bad','2iki':' *Bad','2pdn':'*Good','1pwm':'Bad','3ry2':'Good','3rdo':'Good'}
    
    wild_type = ['1a4h', '1nje', '1syh', '1k3l', '1oh0', '3fs6', '2pdg', '1pwm', '2fzd', '3ry2','3kgu', '3neo', '2iki', '2fr3', '2pzn' ]
    good_aff =['1syi', '3rdo', '3nex', '1ydk', '2g78','3cwk', '3kgt','1v2k', '2pdy', '3f8y', '2pdn', '3mh4', '2pdp']

    c = ['bad_residues','local_residues','ratio','manual', 'affinity']
    #questi cinque dataframe conterranno i risultati di tutte le proteine per ogni threshold scelto
    df_ratio_rmsd075 = pd.DataFrame(index=proteins, columns = c)
    df_ratio_rmsd1 = pd.DataFrame(index=proteins, columns = c)
    df_ratio_rmsd125 = pd.DataFrame(index=proteins, columns = c)
    df_ratio_rmsd15 = pd.DataFrame(index=proteins, columns = c)
    df_ratio_rmsd175 = pd.DataFrame(index=proteins, columns = c)
    
    df_ratio_std025 = pd.DataFrame(index=proteins, columns = c)
    df_ratio_std05 = pd.DataFrame(index=proteins, columns = c)
    df_ratio_std075 = pd.DataFrame(index=proteins, columns = c)
    
    df_ratio_both025 = pd.DataFrame(index=proteins, columns = c)
    df_ratio_both05 = pd.DataFrame(index=proteins, columns = c)
    
    c2=['glob_std', 'glob_std_sc', 'glob_std_bb', 'loc_std', 'loc_std_sc', 'loc_std_bb', 'manual', 'affinity']
    c3=['glob_rmsd', 'glob_rmsd_sc', 'glob_rmsd_bb', 'loc_rmsd', 'loc_rmsd_sc', 'loc_rmsd_bb', 'manual', 'affinity']    

    mis_std = pd.DataFrame(index=proteins, columns = c2)
    mis_rmsd= pd.DataFrame(index=proteins, columns = c3)
    
    #scorro tutte le proteine
    for p_name in proteins:
        print(p_name)
                
        #--------------------------------------Preparazione della proteina--------------------------------------------------------------
        path_af = proteins_folder + '\\' + p_name + '\\' + p_name + '_af.pdb'
        mut = proteins_folder+ '\\' + p_name + '\\mut'
        path_exp = mut + '\\' + os.listdir(mut)[0]
        path_alignment = proteins_folder + '\\' + p_name + '\\' + p_name + '_alignment.pdb'
        #creo il file alignment
        #writeScript(path_exp, path_af, path_alignment, 'script.cxc')
        
        #leggo il file alignment e creo il dataframe af
        ppdb.read_pdb(path_alignment) 
        df_af= ppdb.df['ATOM']
        #leggo il file exp e creo il dataframe exp
        ppdb.read_pdb(path_exp) 
        df_exp_temp1 = ppdb.df['ATOM']

        #Prendo le catene giuste di df_exp
        df_exp_temp2 = df_exp_temp1.loc[df_exp_temp1.chain_id == '']
        if(len(df_exp_temp2)==0):
            df_exp_temp2 = df_exp_temp1.loc[df_exp_temp1.chain_id == 'A']
        if(len(df_exp_temp2)==0):
            df_exp_temp2 =df_exp_temp1.loc[df_exp_temp1.chain_id == 'T']
        
        df_exp = df_exp_temp2.loc[(df_exp_temp2.alt_loc == '')|(df_exp_temp2.alt_loc == 'A')]
        if (len(df_exp)==0):
            df_exp =  df_exp_temp2.loc[df_exp_temp2.alt_loc == 'B']
    
        #elimino gli atomi di idrogeno
        df_af = eliminateH(df_af)
        df_exp = eliminateH(df_exp)
        
        #ordino per numero di residuo e resetto l'indice
        df_af.sort_values(by=['residue_number', 'atom_name'], inplace=True)
        df_af.reset_index(inplace=True, drop=True)
        df_exp.sort_values(by=['residue_number', 'atom_name'], inplace=True)
        df_exp.reset_index(inplace=True, drop=True)
      
        residues_af = list(df_af['residue_number'].drop_duplicates())
        residues_exp = list(df_exp['residue_number'].drop_duplicates())
        if (len(residues_af)!=len(residues_exp)):
            print(len(residues_af)-len(residues_exp))
            for i in range(0, len(residues_exp)):
                print(residues_af[i], residues_exp[i])
                res_af = df_af.loc[df_af.residue_number==residues_af[i]]
                res_exp = df_exp.loc[df_exp.residue_number==residues_exp[i]]
                #if (res_af.iloc[0].loc['residue_name']!=res_exp.iloc[0].loc['residue_name']):
                   # print(res_af)
                    #print(res_exp)
        
        assert(len(residues_af)==len(residues_exp))
        
        diff = len(df_af)-len(df_exp)
        #print(diff)
        for i in range(0, len(residues_exp)):
            res_af = df_af.loc[df_af.residue_number==residues_af[i]]
            res_exp = df_exp.loc[df_exp.residue_number==residues_exp[i]]
            if (len(res_af)!=len(res_exp)):
                #print(residues_af[i], residues_exp[i])
                #print(res_exp[['residue_name', 'residue_number', 'atom_name']])
                #print(res_af[['residue_name', 'residue_number', 'atom_name']])
                df_af.drop(index = res_af.index.tolist(), inplace = True)
                df_exp.drop(index= res_exp.index.tolist(), inplace = True)
        
        diff = len(df_af)-len(df_exp)
        #print(diff)
        #ordino per numero di residuo e resetto l'indice
        df_af.sort_values(by=['residue_number', 'atom_name'], inplace=True)
        df_af.reset_index(inplace=True, drop=True)
        df_exp.sort_values(by=['residue_number', 'atom_name'], inplace=True)
        df_exp.reset_index(inplace=True, drop=True)
      
        #creo la tabella join  
        join = get_join(df_af, df_exp)
        #-----------------------------------------------------------------------------------------------------------------------------------------
        
        
        #--------------------------------------Calcoli su tutta la proteina-----------------------------------------------------------------------
        
        ############################################-----calcoli preliminari-------#####################################################################

        #calcolo il global rmsd
        gl_rmsd = rmsd(df_af, df_exp)
        
        #calcolo il global rmsd sulle catene laterali
        gl_rmsd_SC = rmsd(get_sideChains(df_af),get_sideChains(df_exp))
        
        #calcolo il global rmsd sulla backbone
        gl_rmsd_BB = rmsd(get_backbone(df_af),get_backbone(df_exp))
        
        #calcolo la global std 
        gl_std = join['distance'].std()
        
        #calcolo la global std sulle catene laterali
        join_SC = get_sideChains(join)
        gl_std_SC = join_SC['distance'].std()
        
        #calcolo la global std sulla backbone
        join_BB = get_backbone(join)
        gl_std_BB = join_BB['distance'].std()
        
        ############################################-----calcoli sui residui-------#####################################################################
        
        #calcolo le differenze resiuo per residuo su tutta la proteina
        df_residues_total = calculate_differences(join)

        #----------------------------------------------------------------------------------------------------------------------------------------------        
        
        
        #--------------------------------------Calcoli sulla zona della proteina vicina al ligando------------------------------------------------------------------
        
        ############################################-----calcoli preliminari-------#####################################################################
        #estraggo gli atomi vicini al ligando
        residues_nearLigand = zone(path_alignment, path_exp, dist) 
        
        #prendo gli atomi vicini al ligando nella tabellla join
        join_nearLigand = join.loc[join.residue_number_af.isin(residues_nearLigand)]
        
        df_af_nearLigand = join_nearLigand.loc[: , 'record_name_af':'line_idx_af']
        df_af_nearLigand.rename(columns={"x_coord_af": "x_coord", 'y_coord_af' : 'y_coord', 'z_coord_af':'z_coord'}, inplace=True)
        df_exp_nearLigand = join_nearLigand.loc[: , 'record_name_exp':'line_idx_exp']
        df_exp_nearLigand.rename(columns={"x_coord_exp": "x_coord", 'y_coord_exp' : 'y_coord', 'z_coord_exp':'z_coord'}, inplace=True)
        
        #calcolo il local rmsd
        loc_rmsd = rmsd(df_af_nearLigand, df_exp_nearLigand)
        
        #calcolo il local rmsd sulle catene laterali
        loc_rmsd_SC = rmsd(get_sideChains(df_af_nearLigand),get_sideChains(df_exp_nearLigand))
        
        #calcolo il local rmsd sulla backbone
        loc_rmsd_BB = rmsd(get_backbone(df_af_nearLigand),get_backbone(df_exp_nearLigand))
                
        #calcolo la local std 
        loc_std = join_nearLigand['distance'].std()
        
        #calcolo la local std sulle catene laterali
        join_nearLigand_SC = get_sideChains(join_nearLigand)
        loc_std_SC = join_nearLigand_SC['distance'].std()
        
        #calcolo la local std sulle catene laterali
        join_nearLigand_BB = get_backbone(join_nearLigand)
        loc_std_BB = join_nearLigand_BB['distance'].std()
        
        #calcolo quali sono i residui vicini al ligando nella tabella df_residues_total aggiungendo una colonna di booleani
        residues_nearLigand = [str(n) for n in residues_nearLigand]
        list_bool = df_residues_total.index.isin(residues_nearLigand)
        df_residues_total['IS_NEAR_LIGAND'] = list_bool
        
        ############################################-----calcoli sui residui-------#####################################################################
        
        #calcolo le differenze residuo per residuo solo per i residui vicini al ligando
        df_residues_nearLigand = df_residues_total.loc[list_bool]        
        
        #------------------------------------------------------------------------------------------------------------------------------------------------
        
        #-----------------------------------------------------------dataframes risultati per ogni threshold-------------------------------------------------
        tot = len(df_residues_nearLigand)
        
        affinity = 'No'
        if (p_name in wild_type):
            affinity= '-'
        elif (p_name in good_aff):
            affinity= 'Yes'
        manual = dic[p_name]
        
        #calcolo la ratio a diversi threshold
        flip_RMSD075 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_RMSD>0.75 )]
        f_ratio_RMSD075 = len(flip_RMSD075)/tot
        
        flip_RMSD1 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_RMSD>1 )]
        f_ratio_RMSD1 = len(flip_RMSD1)/tot
        
        flip_RMSD125 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_RMSD>1.25 )]
        f_ratio_RMSD125 = len(flip_RMSD125)/tot
        
        flip_RMSD15 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_RMSD>1.5 )]
        f_ratio_RMSD15 = len(flip_RMSD15)/tot
        
        flip_RMSD175 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_RMSD>1.75 )]
        f_ratio_RMSD175 = len(flip_RMSD175)/tot
        
        
        flip_STD025 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_STD>0.25 )]
        f_ratio_STD025 = len(flip_STD025)/tot
        
        flip_STD05 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_STD>0.5 )]
        f_ratio_STD05 = len(flip_STD05)/tot
        
        flip_STD075 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_STD>0.75 )]
        f_ratio_STD075 = len(flip_STD075)/tot
        
        flip_both025 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_STD>0.25)&(df_residues_nearLigand.SIDECHAIN_RMSD>1.5)]
        f_ratio_both025 = len(flip_both025)/tot
        
        flip_both05 = df_residues_nearLigand.loc[(df_residues_nearLigand.SIDECHAIN_STD>0.5)&(df_residues_nearLigand.SIDECHAIN_RMSD>1.5)]
        f_ratio_both05 = len(flip_both05)/tot
        

        df_ratio_rmsd075.loc[p_name] = [len(flip_RMSD075), tot, round(f_ratio_RMSD075,2), manual, affinity] 
        df_ratio_rmsd1.loc[p_name] = [len(flip_RMSD1), tot, round(f_ratio_RMSD1,2), manual, affinity] 
        df_ratio_rmsd125.loc[p_name] = [len(flip_RMSD125), tot, round(f_ratio_RMSD125,2), manual, affinity] 
        df_ratio_rmsd15.loc[p_name] = [len(flip_RMSD15), tot, round(f_ratio_RMSD15,2), manual, affinity] 
        df_ratio_rmsd175.loc[p_name] = [len(flip_RMSD175), tot, round(f_ratio_RMSD175,2), manual, affinity] 
        
        df_ratio_std025.loc[p_name]=[len(flip_STD025), tot, round(f_ratio_STD025,2), manual, affinity]
        df_ratio_std05.loc[p_name] = [len(flip_STD05), tot, round(f_ratio_STD05,2), manual, affinity] 
        df_ratio_std075.loc[p_name] = [len(flip_STD075), tot, round(f_ratio_STD075,2), manual, affinity] 
        
        df_ratio_both025.loc[p_name] = [len(flip_both025), tot, round(f_ratio_both025,2), manual,affinity] 
        df_ratio_both05.loc[p_name] = [len(flip_both05), tot, round(f_ratio_both05,2), manual,affinity]   
        
        
        mis_std.loc[p_name] = [gl_std, gl_std_SC, gl_std_BB, loc_std, loc_std_SC, loc_std_BB, manual,affinity]
        mis_rmsd.loc[p_name] = [gl_rmsd, gl_rmsd_SC, gl_rmsd_BB, loc_rmsd, loc_rmsd_SC, loc_rmsd_BB, manual, affinity]
        
        #------------------------------------------------------------------------------------------------------------------------------------------------                       
                        
        #--------------------------------------Creo il file di output e scrivo i risultati della singola proteina----------------------------------------
        file_out = open(proteins_folder + '\\' + p_name + '\\' + p_name + '_out.txt', 'w+')
        text = 'Protein name:   ' + p_name + '\n\n'
        
        text += 'Global RMSD:    ' + str(gl_rmsd) + '\n'
        text += 'Global RMSD sidechains:    ' + str(gl_rmsd_SC)+'\n'
        text += 'Global RMSD backbone:    ' + str(gl_rmsd_BB)+'\n\n'
        
        text += 'Local RMSD:    ' + str(loc_rmsd) + '\n'        
        text += 'Local RMSD sidechains: ' + str(loc_rmsd_SC)+'\n'
        text += 'Local RMSD backbone: ' + str(loc_rmsd_BB)+'\n\n'
        
        text += 'Global STD:    ' + str(gl_std) + '\n'
        text += 'Global STD sidechains:    ' + str(gl_std_SC)+'\n'
        text += 'Global STD backbone:    ' + str(gl_std_BB)+'\n\n'
        
        text += 'Local STD:    ' + str(loc_std) + '\n'
        text += 'Local STD sidechains: ' + str(loc_std_SC)+'\n'
        text += 'Local STD backbone: ' + str(loc_std_BB)+'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_RMSD > 0.75\n'
        text += 'Flipped ratio:    ' + str(f_ratio_RMSD075) + '\n'
        text += 'Flipped sidechains:\n' + flip_RMSD075.to_string() +'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_RMSD > 1\n'
        text += 'Flipped ratio:    ' + str(f_ratio_RMSD1) + '\n'
        text += 'Flipped sidechains:\n' + flip_RMSD1.to_string() +'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_RMSD > 1.25\n'
        text += 'Flipped ratio:    ' + str(f_ratio_RMSD125) + '\n'
        text += 'Flipped sidechains:\n' + flip_RMSD125.to_string() +'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_RMSD > 1.5\n'
        text += 'Flipped ratio:    ' + str(f_ratio_RMSD15) + '\n'
        text += 'Flipped sidechains:\n' + flip_RMSD15.to_string() +'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_RMSD > 1.75\n'
        text += 'Flipped ratio:    ' + str(f_ratio_RMSD175) + '\n'
        text += 'Flipped sidechains:\n' + flip_RMSD175.to_string() +'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_STD > 0.25\n'
        text += 'Flipped ratio:    ' + str(f_ratio_STD025) + '\n'
        text += 'Flipped sidechains:\n' + flip_STD025.to_string() +'\n\n\n'

        text += 'Treshold:  SIDECHAIN_STD > 0.5\n'
        text += 'Flipped ratio:    ' + str(f_ratio_STD05) + '\n'
        text += 'Flipped sidechains:\n' + flip_STD05.to_string() +'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_STD > 0.75\n'
        text += 'Flipped ratio:    ' + str(f_ratio_STD075) + '\n'
        text += 'Flipped sidechains:\n' + flip_STD075.to_string() +'\n\n\n'
        
        text += 'Treshold:  SIDECHAIN_STD > 0.25 and SIDECHAIN_RMSD > 1.5\n'
        text += 'Flipped ratio:    ' + str(f_ratio_both025) + '\n'
        text += 'Flipped sidechains:\n' + flip_both025.to_string() +'\n\n\n'  
        
        text += 'Treshold:  SIDECHAIN_STD > 0.5 and SIDECHAIN_RMSD > 1.5\n'
        text += 'Flipped ratio:    ' + str(f_ratio_both05) + '\n'
        text += 'Flipped sidechains:\n' + flip_both05.to_string() +'\n\n\n'   
        
        text += 'Local residues:\n' + df_residues_nearLigand.to_string() + '\n\n'
        text += 'All residues:\n' + df_residues_total.to_string() + '\n\n'
        text += 'All atoms:\n' + join[['residue_number_af', 'residue_name_af', 'atom_name_af', 'residue_name_exp', 'atom_name_exp', 'distance']].to_string()
    
        file_out.write(text)
        file_out.close()
        #------------------------------------------------------------------------------------------------------------------------------------------------
        
    df_ratio_rmsd075.sort_values(by='ratio', inplace=True, ascending=False)
    df_ratio_rmsd1.sort_values(by='ratio', inplace=True, ascending=False)
    df_ratio_rmsd125.sort_values(by='ratio', inplace=True, ascending=False)
    df_ratio_rmsd15.sort_values(by='ratio', inplace=True, ascending=False)
    df_ratio_rmsd175.sort_values(by='ratio', inplace=True, ascending=False)
    
    df_ratio_std025.sort_values(by='ratio', inplace=True, ascending=False)
    df_ratio_std05.sort_values(by='ratio', inplace=True, ascending=False)
    df_ratio_std075.sort_values(by='ratio', inplace=True, ascending=False)
    
    df_ratio_both025.sort_values(by='ratio', inplace=True, ascending=False)
    df_ratio_both05.sort_values(by='ratio', inplace=True, ascending=False)

    results = 'Risultati con threshold rmsd > 0.75\n'
    results+= df_ratio_rmsd075.to_string() + '\n\n\n'
    results+= 'Risultati con threshold rmsd > 1\n'
    results+= df_ratio_rmsd1.to_string() + '\n\n\n'
    results+= 'Risultati con threshold rmsd > 1.25\n'
    results+= df_ratio_rmsd125.to_string() + '\n\n\n'
    results+= 'Risultati con threshold rmsd > 1.5\n'
    results+= df_ratio_rmsd15.to_string() + '\n\n\n'
    results+= 'Risultati con threshold rmsd > 1.75\n'
    results+= df_ratio_rmsd175.to_string() + '\n\n\n'
    
    results+= 'Risultati con threshold std > o.25\n'
    results+= df_ratio_std025.to_string() + '\n\n\n'
    results+= 'Risultati con threshold std > 0.5\n'
    results+= df_ratio_std05.to_string() + '\n\n\n'
    results+= 'Risultati con threshold std > 0.75\n'
    results+= df_ratio_std075.to_string() + '\n\n\n'
        
    results+= 'Risultati con threshold rmsd > 1.5 e std > 0.25\n'
    results+= df_ratio_both025.to_string() + '\n\n\n'
    results+= 'Risultati con threshold rmsd > 1.5 e std > 0.5\n'
    results+= df_ratio_both05.to_string()
    
    file_new = open('file_analizzati.txt', 'w+')
    file_new.write(results)
    file_new.close()
    
    mis_rmsd.sort_values(by='loc_rmsd_sc',ascending=False, inplace=True)
    mis_std.sort_values(by='loc_std_sc',ascending=False, inplace=True)
    file_new = open('tot_protein.txt', 'w+')
    t = 'RMSD su tutta la proteina\n' + mis_rmsd.to_string() + '\n\n'
    t+='\n\n STD su tutta la proteina\n' + mis_std.to_string()
    file_new.write(t)
    file_new.close()
    
    return