"""
Compares output between two TSO500 runs for the same set of samples
"""

import re
from plotnine import *
import pandas as pd
import glob
import os
import scipy
import scipy.stats
import numpy as np
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

def get_sample_names(truth_glob, validation_glob, validation_r_glob):
    truth_samples = glob.glob(truth_glob)
    valid_samples = glob.glob(validation_glob)
    valid_r_samples = glob.glob(validation_r_glob)

    # match all 3 samples by base name and store in df
    # format: base_ID, truth_file, valid_file, valid_file_r
    truth_df = pd.DataFrame(columns=["base_ID", "truth_file"])
    for sample in truth_samples:
        name_components = os.path.basename(sample).split("-")
        base_name = name_components[0] + "-" + name_components[1]
        truth_df.loc[len(truth_df.index)] = [base_name, sample]

    valid_df = pd.DataFrame(columns=["base_ID", "validation_file"])
    for sample in valid_samples:
        name_components = os.path.basename(sample).split("-")
        base_name = name_components[0] + "-" + name_components[1]
        valid_df.loc[len(valid_df.index)] = [base_name, sample]

    valid_r_df = pd.DataFrame(columns=["base_ID", "validation_r_file"])
    for sample in valid_r_samples:
        name_components = os.path.basename(sample).split("-")
        base_name = name_components[0][:-1] + "-" + name_components[1]
        valid_r_df.loc[len(valid_r_df.index)] = [base_name, sample]

    sample_df = pd.merge(left=truth_df, right=valid_df, on="base_ID")
    sample_df = pd.merge(left=sample_df, right=valid_r_df, on="base_ID")

    return sample_df


def make_misc_info_plots(sample_df):
    info_df = pd.DataFrame(columns=["base_ID",
                                    "Total TMB Truth",
                                    "Total TMB Validation",
                                    "Total TMB Validation R",
                                    "Coding Region Size in Megabases Truth",
                                    ("Coding Region Size in "
                                     + "Megabases Validation"),
                                    ("Coding Region Size in "
                                     + "Megabases Validation R"),
                                    ("Number of Passing Eligible "
                                     + "Variants Truth"),
                                    ("Number of Passing Eligible "
                                     + "Variants Validation"),
                                    ("Number of Passing Eligible "
                                     + "Variants Validation R"),
                                    "Usable MSI Sites Truth",
                                    "Usable MSI Sites Validation",
                                    "Usable MSI Sites Validation R",
                                    "Total MSI Sites Unstable Truth",
                                    "Total MSI Sites Unstable Validation",
                                    "Total MSI Sites Unstable Validation R",
                                    "Percent Unstable MSI Sites Truth",
                                    "Percent Unstable MSI Sites Validation",
                                    "Percent Unstable MSI Sites Validation R"])

    # get TMB and MSI values for each sample and sample type

    for index, row in sample_df.iterrows():
        info_truth = get_misc_info(row["truth_file"])
        info_valid = get_misc_info(row["validation_file"])
        info_valid_r = get_misc_info(row["validation_r_file"])
        info_df.loc[len(info_df.index)] = [row["base_ID"],
                                           info_truth[0],
                                           info_valid[0],
                                           info_valid_r[0],
                                           info_truth[1],
                                           info_valid[1],
                                           info_valid_r[1],
                                           info_truth[2],
                                           info_valid[2],
                                           info_valid_r[2],
                                           info_truth[3],
                                           info_valid[3],
                                           info_valid_r[3],
                                           info_truth[4],
                                           info_valid[4],
                                           info_valid_r[4],
                                           info_truth[5],
                                           info_valid[5],
                                           info_valid_r[5]]
        info_df.to_csv("output/TMB_MSI_metrics.tsv", sep="\t")


def make_scatter_plots(sample_df):
    for index, row in sample_df.iterrows():

        venn_name = "output/{}_truth_vs_valid_venn.png".format(row["base_ID"])
        compare_outputs(row["truth_file"],
                        row["validation_file"],
                        "output/{}_intersect.tsv".format(row["base_ID"]),
                        "output/{}_truth_only.tsv".format(row["base_ID"]),
                        "output/{}_valid_only.tsv".format(row["base_ID"]),
                        venn_name)
        truth_name = os.path.basename(row["truth_file"]).split("_")[0]
        valid_name = os.path.basename(row["validation_file"]).split("_")[0]
        valid_r_name = os.path.basename(row["validation_r_file"]).split("_")[0]
        plot_results_scatter("output/{}_intersect.tsv".format(row["base_ID"]),
                             row["base_ID"],
                             "truth_vs_validation",
                             truth_name, valid_name)

        venn_name = ("output/{}R_truth_vs_valid_r_venn.png"
                     .format(row["base_ID"]))
        compare_outputs(row["truth_file"],
                        row["validation_r_file"],
                        "output/{}R_intersect.tsv".format(row["base_ID"]),
                        "output/{}R_truth_only.tsv".format(row["base_ID"]),
                        "output/{}R_valid_only.tsv".format(row["base_ID"]),
                        venn_name)
        plot_results_scatter("output/{}R_intersect.tsv"
                             .format(row["base_ID"]),
                             row["base_ID"],
                             "truth_vs_validation_repeat",
                             truth_name, valid_r_name)


def get_misc_info(file_name):
    file_name = glob.glob(file_name)[0]
    tmb_data = extract_data(file_name, r"^\[TMB\]")
    msi_data = extract_data(file_name, r"^\[MSI\]")

    total_TMB = tmb_data["Total TMB"]
    coding_region = tmb_data["Coding Region Size in Megabases"]
    eligable_variants = tmb_data["Number of Passing Eligible Variants"]
    usable_MSI = msi_data["Usable MSI Sites"]
    MSI_unstable = msi_data["Total MSI Sites Unstable"]
    MSI_percent_unsatable = msi_data["Percent Unstable MSI Sites"]
    return ([total_TMB, coding_region, eligable_variants, usable_MSI,
            MSI_unstable, MSI_percent_unsatable])


def plot_results_venn(unique_truth, unique_valid, intersect, venn_name):
    # plot venn diagram of shared variants and unique to each variant
    venn2(subsets=(len(unique_truth), len(unique_valid), len(intersect)),
          set_labels=('truth', 'val', 'shared'))
    title=venn_name.partition("/")[2].partition("_")[0] + ' Truth vs Validation'
    plt.title(title, fontweight='bold', fontsize=10, pad=25)
    plt.savefig(venn_name)
    plt.close()


def plot_results_scatter(intersect_name, base_ID, comparison_type,
                         truth_ID, valid_ID):
    dat = pd.read_csv(intersect_name, sep='\t')
    # the c.notation was used as a match but this can include intronic
    # regions so remove where p.notation is null
    dat = dat.dropna(subset=['P-Dot Notation_x'])
    dat = dat.dropna(subset=['P-Dot Notation_y'])
    x_lab = "Truth ({})".format(truth_ID)
    y_lab = "Validation ({})".format(valid_ID)

    r = scipy.stats.pearsonr(dat["Allele Frequency_x"],
                             dat["Allele Frequency_y"])[0]
    pearson_r_text = "Pearson R: " + str(np.round(r, 2))

    x_pos = max(dat["Allele Frequency_x"]) * 0.1
    y_pos = max(dat["Allele Frequency_y"]) * 0.9

    plot = (ggplot(dat, aes("Allele Frequency_x", "Allele Frequency_y"))
            + geom_point()
            + stat_smooth(method="lm")
            + geom_text(label=pearson_r_text,
                        x=x_pos, y=y_pos, stat='identity')
            + labs(x=x_lab, y=y_lab)
            + theme_classic()
            + xlim(0, 1)
            + ylim(0, 1))

    output_name = "output/{}_{}_AF.png".format(base_ID, comparison_type)
    plot.save(output_name, height=6, width=10)

    r = scipy.stats.pearsonr(dat["Depth_x"], dat["Depth_y"])[0]
    pearson_r_text = "Pearson R: " + str(np.round(r, 2))

    x_pos = max(dat["Depth_x"]) * 0.1
    y_pos = max(dat["Depth_y"]) * 0.9

    plot = (ggplot(dat, aes("Depth_x", "Depth_y"))
            + geom_point()
            + stat_smooth(method="lm")
            + geom_text(label=pearson_r_text,
                        x=x_pos, y=y_pos, stat='identity')
            + labs(x=x_lab, y=y_lab)
            + theme_classic()
            + xlim(0, max(dat['Depth_x']))
            + ylim(0, max(dat['Depth_y'])))

    output_name = "output/{}_{}_depth.png".format(base_ID, comparison_type)
    plot.save(output_name, height=6, width=10)


def compare_outputs(truth_file, validation_file,
                    output_name,
                    output_name_truth,
                    output_name_valid,
                    venn_name):
    truth_file = glob.glob(truth_file)[0]
    validation_file = glob.glob(validation_file)[0]
    truth_vars = parse_small_variants(truth_file)
    truth_vars.dropna(subset="P-Dot Notation", inplace=True)
    valid_vars = parse_small_variants(validation_file)
    valid_vars.dropna(subset="P-Dot Notation", inplace=True)

    intersect = pd.merge(truth_vars, valid_vars, how="inner",
                         left_on="C-Dot Notation", right_on="C-Dot Notation")
    intersect.to_csv(output_name, sep="\t")

    full = pd.merge(truth_vars, valid_vars,
                    on="C-Dot Notation", how="left", indicator=True)
    df = full.loc[full["_merge"] == "left_only", "C-Dot Notation"]
    truth_only = truth_vars[truth_vars["C-Dot Notation"].isin(df)]
    truth_only.to_csv(output_name_truth, sep="\t")

    full = pd.merge(valid_vars, truth_vars,
                    on="C-Dot Notation", how="left", indicator=True)
    df = full.loc[full["_merge"] == "left_only", "C-Dot Notation"]
    valid_only = valid_vars[valid_vars["C-Dot Notation"].isin(df)]
    valid_only.to_csv(output_name_valid, sep="\t")

    # generate venn diagram
    plot_results_venn(truth_only, valid_only, intersect, venn_name)


def parse_small_variants(file):
    base_name = os.path.basename(file)
    output_name = "temp/" + base_name + ".temp"
    with open(file, "r") as f:
        data_found = False
        with open(output_name, "w") as w:
            for line in f:
                if data_found:
                    if not line.isspace():
                        w.write(line)
                elif re.match(r"^\[Small Variants\]", line):
                    data_found = True
        if not data_found:
            raise Exception("Could not find small variants in {}"
                            .format(file))

    return pd.read_csv(output_name, sep="\t")


def barplots(dat, type):
    df = dat[["base_ID",type+"Truth",
                                        type+"Validation",
                                        type+"Validation R",]]
    df.columns = [col.replace(type, '') for col in df.columns]
    df_long = df.reset_index()
    df_long = pd.melt(df_long,id_vars='base_ID',
                            value_vars=['Truth', 'Validation', 'Validation R'])

    df_long = df_long.rename(columns={'variable': 'Sample type',
                                    'value': type,
                                        'base_ID': 'Sample ID'})


    plot = (ggplot(df_long, aes(x='Sample ID', y=type, fill='Sample type'))
    + geom_col(stat='identity', position='dodge')
    + labs(title=type)
    + theme_classic() 
    + theme(axis_text_x=element_text(rotation=-15, hjust=0.1))
    )

    if type == "Percent Unstable MSI Sites ":
        plot = plot + ylim(0, 100) + labs(title="Percent (%) Unstable MSI Sites ")

    outfile= "output/" + type+" barplots.png"
    plot.save(outfile, height=6, width=10)


def extract_data(file_name, category_regex):
    return_lines = []
    with open(file_name, "r") as f:
        data_found = False
        for line in f:
            if data_found:
                if line.isspace():
                    break
                else:
                    return_lines.append(line)
            elif re.match(category_regex, line):
                data_found = True
        if not data_found:
            raise Exception("Could not find {} in {}"
                            .format(category_regex, file_name))

    return_dict = dict()
    for line in return_lines:
        split = line.split("\t")
        return_dict[split[0]] = split[1]

    return return_dict


if __name__ == "__main__":
    validation_file = (r"dna_cvo_qc_to_compare/"
                       + r"*PCD*"
                       + r"_CombinedVariantOutput.tsv")
    truth_file = (r"dna_cvo_qc_to_compare/"
                  + r"*TSOD*"
                  + r"_CombinedVariantOutput.tsv")
    validation_r_file = (r"dna_cvo_qc_to_compare/"
                         + r"*R-*PC*"
                         + r"_CombinedVariantOutput.tsv")
    sample_df = get_sample_names(truth_file, validation_file,
                                 validation_r_file)
    make_scatter_plots(sample_df)
    make_misc_info_plots(sample_df)

    file_name = 'output/TMB_MSI_metrics.tsv'
    dat = pd.read_csv(file_name, sep='\t')
    barplots(dat, "Total TMB ")
    barplots(dat, "Coding Region Size in Megabases ")
    barplots(dat, "Number of Passing Eligible Variants ")
    barplots(dat, "Usable MSI Sites ")
    barplots(dat, "Total MSI Sites Unstable ")
    barplots(dat, "Percent Unstable MSI Sites ")
    # clean up temp folder
    os.system('rm -rf temp/*')
