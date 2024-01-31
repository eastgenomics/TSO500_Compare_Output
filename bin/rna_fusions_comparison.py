### Tabulate the fusions called from all samples
import pandas as pd
import glob

file_list = glob.glob("rna_fusions_to_compare/*.abridged.tsv.coding_effect")
dfs = []
for file in file_list:
    df = pd.read_csv(file, sep='\t')
    df['sample_name'] = file.split("/")[-1].split(".")[0]
    # shift column 'filename' to first position 
    first_column = df.pop('sample_name') 
    # insert column using insert(position,column_name, 
    # first_column) function 
    df.insert(0, 'sample_name', first_column) 
    dfs.append(df)

# concat the files and save
(pd.concat(dfs,
           ignore_index=True # optional, if you want to keep the index
          )
   .to_csv('output/merged_abdriged_fusions.tsv', sep='\t',
           index=False # optional, if you don't want the index in the output
           )
)

df = pd.concat(dfs, ignore_index=True)


## the clinical samples have the fusions in the CVO file, so get the fusion info from there
file_list = glob.glob("rna_fusions_to_compare/*Output.tsv")
dfs_cvo = []
fusiondict_all = []
for file in file_list:
    # give it random column names so pd can read it in, the names mean nothing
    df_cvo = pd.read_csv(file,names=["V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11"], sep ="\t")
    # figure out rows which contains fusion info
    start_row = df_cvo.index.get_loc(df_cvo.loc[df_cvo['V1'] == '[Fusions]'].index[0])
    end_row= df_cvo.index.get_loc(df_cvo.loc[df_cvo['V1'] == '[Small Variants]'].index[0])
    fusionscvo = list(df_cvo.iloc[int(start_row)+2:int(end_row)-1].V1)
    fusionscvo =  ', '.join(map(str, fusionscvo)) # convert list into comma sep string
    # read into a dictionary
    sample_name = file.split("/")[-1].split("_")[0]
    fusiondict = dict(sample = sample_name, fusions = fusionscvo)
    fusiondict_all.append(fusiondict)

(pd.DataFrame.from_dict(fusiondict_all
          )
   .to_csv('output/merged_cvo_fusions.tsv', sep='\t',
           index=False # optional, if you don't want the index in the output
           )
)
