## this script will pull out the supporting reads for a fusion (=^・^=)
import pandas as pd
from plotnine import *
import collections



def reformat_dataframe(df):
    """  There are many fusion isoforms, we want to report on them
    to distinguish each one we will use the break points as
    the unique identifier along with the fusion name

    Args:
        df : panda dataframe
    """
    # convert the "--" to ":"
    df['#FusionName'] = df['#FusionName'].replace("--", "::", regex=True)
    # make the fusion_id column combining fusiong + breakpoints
    df['fusion_id'] = df['#FusionName'] + "-" + df['LeftBreakpoint'] + "-" +  df['RightBreakpoint']
    df_filt = df[['sample_name', '#FusionName', 'fusion_id', 'JunctionReadCount', 'SpanningFragCount']]
    # filter out clinical samples, they all start with 23 in second column
    df_filt = df_filt[~df_filt['sample_name'].str.contains('-232', regex=False)]
    df_filt = df_filt[~df_filt['sample_name'].str.contains('-233', regex=False)]
    # get the first sample ID
    df_filt['sample_ID'] = df_filt['sample_name'].str.split("-").str[0]
    # get the batch ID info
    df_filt['batch_ID'] = df_filt['sample_name'].str.split("-").str[2]
    # create empty column of what type of dilution or dataset sample is from
    df_filt['type'] = ""
    df_filt['type'] = list(map(lambda x: x.endswith('R'),df_filt["sample_ID"]))
    df_filt['type'] = df_filt['type'].replace(True, '1:4 dilution 23PCR1')
    df_filt.loc[df_filt["batch_ID"] == "24PCR1", "type"] = "1:2 dilution 24PCR1"
    df_filt['type'] = df_filt['type'].replace(False, '1:2 dilution 23PCR1')
    df_filt['sample_ID'] = df_filt['sample_ID'].str.replace(r'R$', '', regex=True)
    # get list of samples
    samples = list(set(df_filt['sample_ID']))

    return df_filt, samples



def barplots(df_filt, samples, metric):
    for sam in samples:
        print(sam)
        df_filt_subset = df_filt[df_filt['sample_ID'].str.contains(sam)]
        df_filt_subset = df_filt_subset.rename(columns={"#FusionName": "FusionName"})
        df_filt_subset['fusion_type'] = df_filt_subset.FusionName + "-" + df_filt_subset.type
        # get list of fusions with type appended and unique
        fusiontype = list(set(df_filt_subset['fusion_type'].tolist()))
        fusiontype_list = [i.split('-', 1)[0] for i in fusiontype]
        # make dictionary of count of fusions
        d = {x:fusiontype_list.count(x) for x in fusiontype_list}
        # filter for fusions that are present in all datasets
        num_of_types = list(set(df_filt_subset['type'].tolist()))
        common_fusions = [k for k, v in d.items() if v >= len(num_of_types)]
        df_filt_subset_common = df_filt_subset.loc[df_filt_subset['FusionName'].isin(common_fusions)]
        base_plot = [
            aes(x='fusion_id', y=metric, fill="type"),
            geom_col(stat='identity', position='dodge'),
            coord_flip(),
            labs(x = "", y=metric, title = sam),
            theme(figure_size=(15, 15))
        ]
        yield ggplot(df_filt_subset_common) + base_plot




def main():
    df = pd.read_csv("output/merged_abdriged_fusions.tsv", sep = "\t")
    df_filt, samples = reformat_dataframe(df)
    print("Creating barplots")
    save_as_pdf_pages(barplots(df_filt, samples, 'JunctionReadCount'),filename = "output/{}_barplots.pdf".format('JunctionReadCount'))
    save_as_pdf_pages(barplots(df_filt, samples, 'SpanningFragCount'),filename = "output/{}_barplots.pdf".format('SpanningFragCount'))



if __name__ == '__main__':
    main()