import glob
import numpy as np
import pandas as pd
from plotnine import *

# Bit of context (◕‿◕✿)
# I am trying to compare the fusions of two samples groups using
# two datasets. The sample groups are: 1) clinical samples and
# 2) validation samples
# There are two files that have fusions: 1) CVO files and
# 2) abridged.tsv.coding_effect files.


def merge_abdriged_files(abriged_path, output_filename):
    """ get all the abriged tsv outputted from pancan pipeline and merge
   into one megatable

    Args:
        abriged_path (str): folder path to where the files are located
    """
    # get all the abriged tsv outputted from pancan pipeline and merge
    # into one megatable
    file_list = glob.glob(abriged_path)
    dfs = []
    for file in file_list:
        df = pd.read_csv(file, sep='\t')
        samplename = file.split("/")[-1].split(".")[0]
        if samplename == "SeraSeqControl":
            samplename = "SSC"
        elif samplename.count("-") == 1:
            samplename = samplename.replace("-", "")
        df['sample_name'] = samplename
        first_column = df.pop('sample_name')
        df.insert(0, 'sample_name', first_column)
        dfs.append(df)

    # concat the files and save
    (pd.concat(dfs, ignore_index=True)
        .to_csv('output/{}'.format(output_filename), sep='\t', index=False))


def merge_cvo_files(cvo_path, output_filename):
    """ get all the cvo tsv outputted from tso500 pipeline and merge
   into one megatable

    Args:
        abriged_path (str): folder path to where the files are located
    """
    file_list = glob.glob(cvo_path)
    dfs_cvo = []
    fusiondict_all = []
    for file in file_list:
        # give it random column names so pd can read it in,
        # the names mean nothing
        df_cvo = pd.read_csv(file, names=["V1", "V2", "V3", "V4", "V5",
                                          "V6", "V7", "V8", "V9", "V10",
                                          "V11"], sep="\t")
        # figure out rows which contains fusion info
        start_row = df_cvo.index.get_loc(
            df_cvo.loc[df_cvo['V1'] == '[Fusions]'].index[0]
            )
        end_row = df_cvo.index.get_loc(
            df_cvo.loc[df_cvo['V1'] == '[Small Variants]'].index[0]
            )
        fusionscvo = list(df_cvo.iloc[int(start_row)+2:int(end_row)-1].V1)
        # convert list into comma sep string
        fusionscvo = ', '.join(map(str, fusionscvo))
        # read into a dictionary
        sample_name = file.split("/")[-1].split("_")[0]
        fusiondict = dict(sample=sample_name, fusions=fusionscvo)
        fusiondict_all.append(fusiondict)

    (pd.DataFrame.from_dict(fusiondict_all).to_csv(
        'output/{}'.format(output_filename), sep='\t', index=False))


def find_clinical_sample_fusions(merged_abridged,
                                  merged_cvo,
                                  output_filename):
    """For the clinical samples (second field has 23) merge all the
    fusions called from the merged_cvo & merged_abriged.tsv reports.

    First make dictionaries for fusions from the merged abridged and
    then from merged cvo. Secondly make a dataframe from both
    dictionaries and save it.

    Args:
        merged_abridged (pandas df): merged abriged fusions file
        generated from function merge_abdriged_files
        merged_cvo (pandas df): merged cvo fusions file
        generated from function merge_cvo_files
        output_filename (str): output file name
    """
    # 1. Make a dictionary of the fusions
    # from pancan pipeline merged_abriged.tsv
    merged_ab = pd.read_csv('output/{}'.format(merged_abridged), sep='\t')
    merged_ab = merged_ab[["sample_name", "#FusionName"]]
    # filter for clinical samples from abridged df, they all start
    # with 23 in second column
    merged_ab_cs = merged_ab.set_index(
                                    'sample_name'
                                    ).filter(
                                    regex='-23', axis=0
                                    )
    c_samples = list(set(merged_ab_cs.index))
    merged_ab_cs['sample_name'] = merged_ab_cs.index
    pancan_fusions = []
    for sam in c_samples:
        temp_dict = {}
        temp_dict['sample'] = sam
        sam_df = merged_ab[merged_ab['sample_name'] == sam]
        fusions = list(set(sam_df['#FusionName'].tolist()))
        temp_dict['fusions'] = fusions
        pancan_fusions.append(temp_dict)

    # 2. Make a dictionary of the fusions from the TSO500 pipeline merged_cvo
    merged_cvo = pd.read_csv('output/{}'.format(merged_cvo), sep='\t')
    # filter for clinical samples from abridged df,
    # they all start with 23 in second column
    merged_cvo_cs = merged_cvo.set_index('sample').filter(regex='-23', axis=0)
    c_samples = list(set(merged_cvo_cs.index))
    merged_cvo_cs['sample'] = merged_cvo_cs.index
    tso_fusions = []
    for sam in c_samples:
        temp_dict = {}
        temp_dict['sample'] = sam
        sam_df = merged_cvo_cs[merged_cvo_cs['sample'] == sam]
        fusions = list(set(sam_df['fusions'].tolist()))
        temp_dict['fusions'] = fusions
        tso_fusions.append(temp_dict)

    # convert dictionaries into dataframes and add a column to show which
    # pipeline they were processed through and convert list to string
    pancan_fusions_df = pd.DataFrame.from_dict(pancan_fusions)
    pancan_fusions_df['fusions'] = (
        pancan_fusions_df['fusions']
        .transform(
            lambda x: ", ".join(map(str, x))
        )
    )
    tso_fusions_df = pd.DataFrame.from_dict(tso_fusions)
    tso_fusions_df['fusions'] = (
        tso_fusions_df['fusions']
        .transform(
            lambda x: ", ".join(map(str, x))
        )
    )

    df = pd.concat([pancan_fusions_df, tso_fusions_df], ignore_index=True)
    df = df.sort_values(by=['sample'])

    # create final table with the clinical samples and
    # fusions called from both pipelines
    df = pd.DataFrame(columns=['sample', 'pancan', 'tso500'])
    df['sample'] = c_samples
    # for each samples, get the fusions from the pancan dataframe
    # and tso500 dataframe. If empty or "nan", then skip
    for i in range(len(df)):
        sam = df.loc[i, 'sample']
        if sam in list(pancan_fusions_df['sample']):
            val = pancan_fusions_df.loc[
                pancan_fusions_df['sample'] == sam
                ].fusions
            if val is not val.empty:
                sam_fusion = val.values[0]
                df.loc[i, 'pancan'] = sam_fusion
        if "23TSO" in sam:
            df.loc[i, 'pancan'] = ""
        if sam in list(tso_fusions_df['sample']):
            val = tso_fusions_df.loc[tso_fusions_df['sample'] == sam].fusions
            if val is not val.empty:
                if list(val)[0] != "nan":
                    sam_fusion = val.values[0]
                    df.loc[i, 'tso500'] = sam_fusion

    df.to_csv('output/{}'.format(output_filename), sep="\t")


def find_validation_sample_fusions(merged_abridged,
                                  output_fusions):
    """For the validation samples (missing 23 in the second field)
    merge all the fusions called from merged_abriged.tsv reports.

    First make dictionaries for fusions from the merged abridged and
    then make a dataframe from the dictionary. Then categorise the
    samples based on the experiment type.

    Args:
        merged_abridged (pandas df): merged abriged fusions file
        generated from function merge_abdriged_files
        output_filename (str): output file name
    """
    # Now lets focus on validation samples.
    # We can read all in and merge the fusions
    merged_ab = pd.read_csv('output/{}'.format(merged_abridged), sep='\t')
    merged_ab = merged_ab[["sample_name", "#FusionName"]]
    # filter for clinical samples from abridged df,
    # they all start with 23 in second column. If we search for just
    # 23 it can filter out the batch 23PCR so we want just the
    # samples with second field starting 232 or 233
    merged_ab = merged_ab[
        ~merged_ab['sample_name'].str.contains('-232', regex=False)
        ]
    merged_ab = merged_ab[
        ~merged_ab['sample_name'].str.contains('-233', regex=False)
        ]
    c_samples = list(set(merged_ab['sample_name']))
    pancan_fusions = []
    for sam in c_samples:
        temp_dict = {}
        temp_dict['sample'] = sam
        sam_df = merged_ab[merged_ab['sample_name'] == sam]
        fusions = list(set(sam_df['#FusionName'].tolist()))
        temp_dict['fusions'] = fusions
        pancan_fusions.append(temp_dict)

    # convert dictionaries into dataframes and add a column to show which
    # pipeline they were processed through and convert list to string
    pancan_fusions_df = pd.DataFrame.from_dict(pancan_fusions)
    pancan_fusions_df['fusions'] = (
        pancan_fusions_df['fusions']
        .transform(
            lambda x: ", ".join(map(str, x))
        )
    )

    pancan_fusions_df = pancan_fusions_df.sort_values(by=['sample'])
    # get the first sample ID
    pancan_fusions_df['sample_ID'] = pancan_fusions_df[
        'sample'
        ].str.split("-").str[0]
    # get the batch ID info
    pancan_fusions_df['batch_ID'] = pancan_fusions_df[
        'sample'
        ].str.split("-").str[2]
    # create empty column of what type of dilution or dataset sample is from
    pancan_fusions_df['type'] = ""
    pancan_fusions_df['type'] = list(
        map(
        lambda x: x.endswith('R'), pancan_fusions_df["sample_ID"]
        )
        )
    pancan_fusions_df['type'] = pancan_fusions_df['type'].replace(
        True, '1:4 dilution 23PCR1'
        )
    pancan_fusions_df.loc[pancan_fusions_df["batch_ID"] == "24PCR1", "type"] = "1:2 dilution 24PCR1"
    pancan_fusions_df['type'] = pancan_fusions_df['type'].replace(
        False, '1:2 dilution 23PCR1'
        )
    # if sample name is just one, then its the neat samples
    for i in range(len(pancan_fusions_df)):
        count = pancan_fusions_df.loc[i, 'sample'].count("-")
        if count <= 1:
            pancan_fusions_df.loc[i, 'type'] = "Neat 22"
            pancan_fusions_df.loc[i, 'batch_ID'] = "Neat 22"
    pancan_fusions_df['sample_ID'] = pancan_fusions_df[
        'sample_ID'
        ].str.replace(
            r'R$', '', regex=True
            )
    pancan_fusions_df = pancan_fusions_df[[
        "sample", "sample_ID", "batch_ID", "type", "fusions"
        ]]
    pancan_fusions_df.to_csv(
        "output/{}".format(output_fusions), sep="\t"
        )

    return pancan_fusions_df


def performance_calculations(pancan_fusions_df, expected_fusions,
                            output_fusions_filename,
                            output_fusions_count_filename):
    """_summary_

    Args:
        pancan_fusions_df (pandas df): df contains fusions of each
        sample generated from function find_validation_sample_fusions
        expected_fusions (pandas df): df where there are
        two columns detailing sample and fusions. Each row being a
        sample with fusions seperated by ", "
        output_fusions_filename (str): output file name
        output_fusions_count_filename (str): output counts file name

    Returns:
        pandas dataframe: dataframe where each row is a sample and metrics
        such as true positives (TP), false positives (FP) and false
        negatives (FN) is listed
    """
    # Now we can read in the expected fusions and match based on name
    expected_fusions = pd.read_csv(expected_fusions,
                                    names=[
                                        "sample_ID", "fusions"
                                        ], sep="\t")

    df = pancan_fusions_df.merge(
                            expected_fusions, on='sample_ID', how='left'
                            )
    df = df.rename(columns={"fusions_x": "observed_fusion", "fusions_y": "expected_fusion"})
    # the observed fusions have two comma seperated so we will repalce
    df['observed_fusion'] = df['observed_fusion'].replace(
        "--", "::", regex=True
        )
    df['expected_fusion'] = df['expected_fusion'].replace(
        "-", "::", regex=True
        )

    # make extra column
    df['TP'] = ""
    df['FP'] = ""
    df['FN'] = ""
    # have another df just to tabulate
    df_tabulate = pd.DataFrame(columns=["sample", "sample_ID", "batch_ID",
                                        "type", "observed_fusion",
                                        "expected_fusion",
                                        "TP", "FP", "FN"])
    df_tabulate["sample"] = df["sample"]
    df_tabulate["sample_ID"] = df["sample_ID"]
    df_tabulate["batch_ID"] = df["batch_ID"]
    df_tabulate["type"] = df["type"]
    df_tabulate["observed_fusion"] = df["observed_fusion"]
    df_tabulate["expected_fusion"] = df["expected_fusion"]

    # we want to populate the TP, FP and FN negatives. For each row,
    # we take take the observed fusions and expected fusions and
    # calculate the TP, FP and FN and put the list in the df and the
    # total count in the df_tabulate object.

    for i in range(len(df)):
        # convert comma sep to lists
        obs_fusion = df.loc[i, 'observed_fusion']
        obs_fusion = obs_fusion.split(", ")
        df_tabulate.loc[i, 'observed_fusion'] = len(obs_fusion)
        exp_fusion = df.loc[i, 'expected_fusion']
        if isinstance(exp_fusion, float):
            exp_fusion = []
            df_tabulate.loc[i, 'expected_fusion'] = 0
        else:
            exp_fusion = exp_fusion.split(", ")
            df_tabulate.loc[i, 'expected_fusion'] = len(exp_fusion)
        # Intersect/Shared
        true_pos = list(set(obs_fusion).intersection(set(exp_fusion)))
        true_pos = [sub.replace(' ', ', ') for sub in true_pos]
        # convert from list to comma seperated string
        true_pos_list = ', '.join(map(str, true_pos)).split(", ")
        # empty fusion is just '', so need to set it as empty list
        if '' in true_pos_list:
            true_pos_list = []
        true_pos = ', '.join(map(str, true_pos))
        df.loc[i, 'TP'] = true_pos
        df_tabulate.loc[i, 'TP'] = len(true_pos_list)
        # only seen in expected fusions
        false_pos = list(set(obs_fusion) - set(exp_fusion))
        false_pos = [sub.replace(' ', ', ') for sub in false_pos]
        # convert from list to comma seperated string
        false_pos_list = ', '.join(map(str, false_pos)).split(", ")
        # empty fusion is just '', so need to set it as empty list
        if '' in false_pos_list:
            false_pos_list = []
        false_pos = ', '.join(map(str, false_pos))
        df.loc[i, 'FP'] = false_pos
        df_tabulate.loc[i, 'FP'] = len(false_pos_list)
        # only seen in observed fusions
        false_neg = list(set(exp_fusion) - set(obs_fusion))
        false_neg = [sub.replace(' ', ', ') for sub in false_neg]
        # convert from list to comma seperated string
        false_neg_list = ', '.join(map(str, false_neg)).split(", ")
        # empty fusion is just '', so need to set it as empty list
        if '' in false_neg_list:
            false_neg_list = []
        false_neg = ', '.join(map(str, false_neg))
        df.loc[i, 'FN'] = false_neg
        df_tabulate.loc[i, 'FN'] = len(false_neg_list)

    # rename headers
    df = df.rename(columns={"sample": "Sample",
                            "sample_ID": "Sample ID",
                            "batch_ID": "Batch_ID",
                            "type": "Type",
                            "observed_fusion": "Observed Fusions",
                            "expected_fusion": "Expected Fusions",
                            "TP": "True Positives",
                            "FP": "False Positives",
                            "FN": "False Negatives"})

    df_tabulate = df_tabulate.rename(columns={
                                        "sample": "Sample",
                                        "sample_ID": "Sample ID",
                                        "batch_ID": "Batch_ID",
                                        "type": "Type",
                                        "observed_fusion": "Observed Fusions",
                                        "expected_fusion": "Expected Fusions",
                                        "TP": "True Positives",
                                        "FP": "False Positives",
                                        "FN": "False Negatives"})

    # all fusions listed here
    df.to_csv(
        "output/{}".format(output_fusions_filename), sep="\t"
        )
    # tabulate with just numbers here and output it for plots later
    df_tabulate.to_csv(
        "output/{}".format(output_fusions_count_filename), sep="\t"
        )

    return df_tabulate


def barplot_performance(fusions_count_df):
    """ The count fusion dataframe can be plotted as barplots with
        different types as groups.
    """
    plot_cols = ["Observed Fusions", "True Positives",
                  "False Positives", "False Negatives"]

    for col in plot_cols:
        print(col)
        fusions_count_df[col] = pd.to_numeric(fusions_count_df[col])
        max_y_axis = int(fusions_count_df[col].max())
        plot = (ggplot(fusions_count_df, aes(
                    x="Sample ID", y=col, fill="Type")
                    )
                    + geom_col(stat='identity', position='dodge')
                    + labs(x="Sample", y=col, title=col)
                    + theme_classic()
                    + theme(
                        axis_text_x=element_text(rotation=-45, hjust=0.1)
                        )
                )

        plot = plot + scale_y_continuous(limits=(0, max_y_axis))
        output_name = "output/fusions_barplots_{}.png".format(col)
        plot.save(output_name, height=6, width=10)


def main():

    merge_abdriged_files(
        abriged_path="rna_fusions_to_compare/*.abridged.tsv.coding_effect",
        output_filename="merged_abdriged_fusions.tsv")

    merge_cvo_files(
        cvo_path="rna_fusions_to_compare/*Output.tsv",
        output_filename="merged_cvo_fusions.tsv")

    find_clinical_sample_fusions(
        merged_cvo="merged_cvo_fusions.tsv",
        merged_abridged="merged_abdriged_fusions.tsv",
        output_filename="clinical_samples_fusions.tsv")

    pancan_fusions_dat = find_validation_sample_fusions(
        merged_abridged="merged_abdriged_fusions.tsv",
        output_fusions="validation_samples_fusions.tsv")

    fusions_counts = performance_calculations(
        pancan_fusions_df=pancan_fusions_dat,
        expected_fusions="rna_fusions_to_compare/expected_fusions.txt",
        output_fusions_filename="validation_samples_and_expected_fusions.tsv",
        output_fusions_count_filename="validation_samples_and_expected_fusions_count.tsv")

    barplot_performance(
        fusions_count_df=fusions_counts
    )


if __name__ == '__main__':
    main()
