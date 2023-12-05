from biopandas.pdb import PandasPdb
import os
import preparazione as pr
import totmisurations as mis
cwd = os.getcwd()
ppdb = PandasPdb()

#La funzione apply applica il metodo, scrivendo i riultati sul file identificato dal parametro "path_out".

def apply(p_name, path_af, path_exp, path_alignment, path_out, dist = 5):

    #--------------------------------------Preparazione della proteina--------------------------------------------------------------
    #creo il file alignment
    pr.writeScript(path_exp, path_af, path_alignment, 'script.cxc')
        
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
    if(len(df_exp_temp2)==0):
        df_exp_temp2 =df_exp_temp1.loc[df_exp_temp1.chain_id == 'X']
        
    df_exp = df_exp_temp2.loc[(df_exp_temp2.alt_loc == '')|(df_exp_temp2.alt_loc == 'A')]
    if (len(df_exp)==0):
        df_exp =  df_exp_temp2.loc[df_exp_temp2.alt_loc == 'B']
    
    #elimino gli atomi di idrogeno
    df_af = pr.eliminateH(df_af)
    df_exp = pr.eliminateH(df_exp)
        
    #ordino per numero di residuo e resetto l'indice
    df_af.sort_values(by=['residue_number', 'atom_name'], inplace=True)
    df_af.reset_index(inplace=True, drop=True)
    df_exp.sort_values(by=['residue_number', 'atom_name'], inplace=True)
    df_exp.reset_index(inplace=True, drop=True)

    residues_af = list(df_af['residue_number'].drop_duplicates())
    residues_exp = list(df_exp['residue_number'].drop_duplicates())

    assert(len(residues_af)==len(residues_exp))
        
    diff = len(df_af)-len(df_exp)
    for i in range(0, len(residues_exp)):
        res_af = df_af.loc[df_af.residue_number==residues_af[i]]
        res_exp = df_exp.loc[df_exp.residue_number==residues_exp[i]]
        if (len(res_af)!=len(res_exp)):
            df_af.drop(index = res_af.index.tolist(), inplace = True)
            df_exp.drop(index= res_exp.index.tolist(), inplace = True)
        
    #ordino per numero di residuo e resetto l'indice
    df_af.sort_values(by=['residue_number', 'atom_name'], inplace=True)
    df_af.reset_index(inplace=True, drop=True)
    df_exp.sort_values(by=['residue_number', 'atom_name'], inplace=True)
    df_exp.reset_index(inplace=True, drop=True)
      
    #creo la tabella join  
    join = pr.get_join(df_af, df_exp)
    #-----------------------------------------------------------------------------------------------------------------------------------------


    #--------------------------------------Calcoli su tutta la proteina-----------------------------------------------------------------------
    #calcoli su tutta la struttura:
    #calcolo il global rmsd
    gl_rmsd = mis.rmsd(df_af, df_exp)
        
     #calcolo il global rmsd sulle catene laterali
    gl_rmsd_SC = mis.rmsd(mis.get_sideChains(df_af), mis.get_sideChains(df_exp))
        
    #calcolo il global rmsd sulla backbone
    gl_rmsd_BB = mis.rmsd(mis.get_backbone(df_af), mis.get_backbone(df_exp))
        
    #calcolo la global std 
    gl_std = join['distance'].std()
        
    #calcolo la global std sulle catene laterali
    join_SC = mis.get_sideChains(join)
    gl_std_SC = join_SC['distance'].std()
        
    #calcolo la global std sulla backbone
    join_BB = mis.get_backbone(join)
    gl_std_BB = join_BB['distance'].std()
        
    #calcoli per ogni residuo:    
    #calcolo le differenze resiuo per residuo su tutta la proteina
    df_residues_total = mis.calculate_differences(join)

    #----------------------------------------------------------------------------------------------------------------------------------------------        
        
    
    #--------------------------------------Calcoli sulla zona della proteina vicina al ligando------------------------------------------------------------------
    #calcoli su tutta la struttura:
    #estraggo gli atomi vicini al ligando
    residues_nearLigand = mis.zone(path_alignment, path_exp, dist) 
        
    #prendo gli atomi vicini al ligando nella tabellla join
    join_nearLigand = join.loc[join.residue_number_af.isin(residues_nearLigand)]
        
    df_af_nearLigand = join_nearLigand.loc[: , 'record_name_af':'line_idx_af']
    df_af_nearLigand.rename(columns={"x_coord_af": "x_coord", 'y_coord_af' : 'y_coord', 'z_coord_af':'z_coord'}, inplace=True)
    df_exp_nearLigand = join_nearLigand.loc[: , 'record_name_exp':'line_idx_exp']
    df_exp_nearLigand.rename(columns={"x_coord_exp": "x_coord", 'y_coord_exp' : 'y_coord', 'z_coord_exp':'z_coord'}, inplace=True)
        
    #calcolo il local rmsd
    loc_rmsd = mis.rmsd(df_af_nearLigand, df_exp_nearLigand)
        
    #calcolo il local rmsd sulle catene laterali
    loc_rmsd_SC = mis.rmsd(mis.get_sideChains(df_af_nearLigand), mis.get_sideChains(df_exp_nearLigand))
        
    #calcolo il local rmsd sulla backbone
    loc_rmsd_BB = mis.rmsd(mis.get_backbone(df_af_nearLigand), mis.get_backbone(df_exp_nearLigand))
                
    #calcolo la local std 
    loc_std = join_nearLigand['distance'].std()
        
    #calcolo la local std sulle catene laterali
    join_nearLigand_SC = mis.get_sideChains(join_nearLigand)
    loc_std_SC = join_nearLigand_SC['distance'].std()
        
    #calcolo la local std sulle catene laterali
    join_nearLigand_BB = mis.get_backbone(join_nearLigand)
    loc_std_BB = join_nearLigand_BB['distance'].std()
        
    #calcolo quali sono i residui vicini al ligando nella tabella df_residues_total aggiungendo una colonna di booleani
    residues_nearLigand = [str(n) for n in residues_nearLigand]
    list_bool = df_residues_total.index.isin(residues_nearLigand)
    df_residues_total['IS_NEAR_LIGAND'] = list_bool
        
    #calcoli per ogni residuo:   
    #calcolo le differenze residuo per residuo solo per i residui vicini al ligando
    df_residues_nearLigand = df_residues_total.loc[list_bool]        
    
    
    #--------------------------------------Creo il file di output e scrivo i risultati della singola proteina----------------------------------------
    file_out = open(path_out, 'w+')
    text = 'AlphaFold2 file:   ' + path_af + '\n'
    text += 'Protein name:  ' + p_name + '\n\n'
        
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
        
    text += 'Local residues:\n' + df_residues_nearLigand.to_string() + '\n\n'
    text += 'All residues:\n' + df_residues_total.to_string() + '\n\n'
    text += 'All atoms:\n' + join[['residue_number_af', 'residue_name_af', 'atom_name_af', 'residue_name_exp', 'atom_name_exp', 'distance']].to_string()
    
    file_out.write(text)
    file_out.close()
    #------------------------------------------------------------------------------------------------------------------------------------------------
    return
