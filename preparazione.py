#Questo file contiene le funzioni che servono per il primo step del metodo, ovvero la preparazione dei dataframe.

import pandas as pd
import os
import subprocess as sp
import numpy as np

#eliminateH(df) prende un dataframe che rappresenta una proteina e restituisce il dataframe senza gli atomi di idrgeno
def eliminateH(df):
    #le righe di temp contengono true solo se l'atomo è un idrogeno
    temp = ~df['atom_name'].str.startswith("H")
    ris = df.loc[temp]
    return ris.reset_index(drop=True)

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
