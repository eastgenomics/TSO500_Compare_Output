### Tabulate the fusions called from all samples
from heapq import merge
from unittest import skip
import pandas as pd
import glob


## get all the abriged tsv outputted from pancan pipeline and merge
## into one megatable
file_list = glob.glob("rna_fusions_to_compare/*.abridged.tsv.coding_effect")
dfs = []
for file in file_list:
    df = pd.read_csv(file, sep='\t')
    df['sample_name'] = file.split("/")[-1].split(".")[0]
    first_column = df.pop('sample_name') 
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
## merge into a megatable
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

## For the clinical samples, lets compare the merged_cvo & merged_abriged.tsv reports
## 1. Make a dictionary of the fusions from pancan pipeline merged_abriged.tsv
merged_ab = pd.read_csv('output/merged_abdriged_fusions.tsv', sep='\t')
merged_ab = merged_ab[["sample_name", "#FusionName"]]
# filter for clinical samples from abridged df, they all start with 23 in second column
merged_ab_cs = merged_ab.set_index('sample_name').filter(regex='-23', axis=0)
c_samples = list(set(merged_ab_cs.index))
merged_ab_cs['sample_name'] = merged_ab_cs.index
pancan_fusions = []
for sam in c_samples:
    temp_dict = {}
    temp_dict['sample'] = sam
    sam_df = merged_ab_cs[merged_ab_cs['sample_name'].str.contains(sam)]
    fusions = list(set(sam_df['#FusionName'].tolist()))
    temp_dict['fusions'] = fusions
    pancan_fusions.append(temp_dict)

## 2. Make a dictionary of the fusions from the TSO500 pipeline merged_cvo
merged_cvo = pd.read_csv('output/merged_cvo_fusions.tsv', sep='\t')
# filter for clinical samples from abridged df, they all start with 23 in second column
merged_cvo_cs = merged_cvo.set_index('sample').filter(regex='-23', axis=0)
c_samples = list(set(merged_cvo_cs.index))
merged_cvo_cs['sample'] = merged_cvo_cs.index
tso_fusions = []
for sam in c_samples:
    temp_dict = {}
    temp_dict['sample'] = sam
    sam_df = merged_cvo_cs[merged_cvo_cs['sample'].str.contains(sam)]
    fusions = list(set(sam_df['fusions'].tolist()))
    temp_dict['fusions'] = fusions
    tso_fusions.append(temp_dict)


# convert dictionaries into dataframes and add a column to show which
# pipeline they were processed through and convert list to string
pancan_fusions_df = pd.DataFrame.from_dict(pancan_fusions)
#pancan_fusions_df['pipeline'] =  pd.Series(["pancan" for x in range(len(pancan_fusions_df.index))])
pancan_fusions_df['fusions'] = (
    pancan_fusions_df['fusions']
    .transform(
        lambda x: ", ".join(map(str,x))
    )
)
tso_fusions_df = pd.DataFrame.from_dict(tso_fusions)
#tso_fusions_df['pipeline'] =  pd.Series(["tso500" for x in range(len(tso_fusions_df.index))])
tso_fusions_df['fusions'] = (
    tso_fusions_df['fusions']
    .transform(
        lambda x: ", ".join(map(str,x))
    )
)

df = pd.concat([pancan_fusions_df,tso_fusions_df], ignore_index=True)
df = df.sort_values(by=['sample'])

if list(pancan_fusions_df['sample']) == list(tso_fusions_df['sample']):
    print('Same rownames')

## create final table with the clinical samples and 
# fusions called from both pipelines
df = pd.DataFrame(columns=['sample','pancan','tso500'])
df['sample'] = c_samples
# for eachh samples, get the fusions from the pancan dataframe
# and tso500 dataframe. If empty or "nan", then skip
for i in range(len(df)):
    sam = df.loc[i,'sample']
    if sam in list(pancan_fusions_df['sample']):
        val = pancan_fusions_df.loc[pancan_fusions_df['sample'] == sam].fusions
        if val is not val.empty:
            sam_fusion = val.values[0]
            df.loc[i,'pancan'] = sam_fusion
    if "23TSO" in sam:
        df.loc[i,'pancan'] = ""
    if sam in list(tso_fusions_df['sample']):
        val = tso_fusions_df.loc[tso_fusions_df['sample'] == sam].fusions
        if val is not val.empty:
            if list(val)[0] != "nan":
                sam_fusion = val.values[0]
                df.loc[i,'tso500'] = sam_fusion

df.to_csv("output/clinical_samples_fusions.tsv", sep="\t")