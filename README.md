# Metodo_per_analisi_di_proteine
Main.py prende in input una cartella e il parametro dist (5 di defaut). 
La cartella deve contenere una sottocartella per ogni proteina che si vuole analizzare, che abbia il nome della proteina. 
Nella cartella di ogni proteina deve esserci il file di AlphaFold2 con il nome che rispetti la struttura: '*nomeproteina*_af.pdb'.
La cartella deve anche contenere una cartella chiamata wt contenente il file wt sperimentale.
Per ogni proteina viene creato un file '_out.txt' contenente i risultati per quela proteina.
