import glob
import numpy as np
import pandas as pd
from plotnine import *

general_stats="rna_qc_to_compare/multiqc_general_stats_*"
output_filename="multiqc_general_stats_all.tsv"
file_list = glob.glob(general_stats)
dfs = []
for file in file_list:
    df = pd.read_csv(file, sep='\t')
    dfs.append(df)

# concat the files and save
(pd.concat(dfs, ignore_index=True)
.to_csv('output/{}'.format(output_filename), sep='\t', index=False))

concant_dfs = pd.concat(dfs, ignore_index=True)
concant_dfs['sample_name'] = concant_dfs["Sample"].str.split("_", expand=True)[0]
concant_dfs['sample_ID'] = concant_dfs["Sample"].str.split("-", expand=True)[0]
concant_dfs['batch_ID'] = concant_dfs["Sample"].str.split("-", expand=True)[2]
## add experimental type
concant_dfs['type'] = ""
concant_dfs['type'] = list(map(lambda x: x.endswith('R'), concant_dfs["sample_ID"]))
concant_dfs['type'] = concant_dfs['type'].replace(True, '1:4 dilution 23PCR1')
concant_dfs.loc[concant_dfs["batch_ID"] == "24PCR1", "type"] = "1:2 dilution 24PCR1"
concant_dfs['type'] = concant_dfs['type'].replace(False, '1:2 dilution 23PCR1')
concant_dfs.loc[concant_dfs["batch_ID"] == "24PCD1", "type"] = "24PCD1"



## rename headers to something clearer
concant_dfs = concant_dfs.rename(columns={
                            'Picard_mqc-generalstats-picard-PCT_MRNA_BASES' : "mRNA_bases",
                            'Picard_mqc-generalstats-picard-PCT_RIBOSOMAL_BASES': "Ribosomal_bases",
                            'RNA-SeQC_mqc-generalstats-rna_seqc-Expression_Profiling_Efficiency': "Expression_Profiling_Efficiency",
                            'RNA-SeQC_mqc-generalstats-rna_seqc-Genes_Detected': "Genes_Detected",
                            'RNA-SeQC_mqc-generalstats-rna_seqc-rRNA_rate': "rRNA_rate",
                            'FastQC_mqc-generalstats-fastqc-total_sequences': "Total_Sequences",
                            'FastQC_mqc-generalstats-fastqc-avg_sequence_length': "Length",
                            'FastQC_mqc-generalstats-fastqc-percent_duplicates': "Duplicates",
                            'FastQC_mqc-generalstats-fastqc-percent_gc': "GC",
                            'FastQC_mqc-generalstats-fastqc-percent_fails': "Fails",})


## plots the fastqc data
df_fastq = concant_dfs[~concant_dfs['Sample'].str.contains(".star")]
cols = ["Total_Sequences", "Length", "Duplicates", "GC", "Fails"]
for col in cols:
    plot = (ggplot(df_fastq, aes(
                x="sample_name", y=col, fill = "type")
                )
                + geom_col(stat='identity', position='dodge')
                + labs(x="",y=col, title=col)
                + coord_flip()
                + theme_classic()
                + theme(
                    axis_text_x=element_text(rotation=-15, hjust=0.1)
                    )
            )
    output_name = "output/qc_fusions_barplots_{}.png".format(col)
    plot.save(output_name, height=10, width=10)

## plots the rnaseq data
df_rnaseqc = concant_dfs[concant_dfs['Sample'].str.contains(".star")]
cols = ["Expression_Profiling_Efficiency", "Genes_Detected", "rRNA_rate", "mRNA_bases"]
for col in cols:
    plot = (ggplot(df_rnaseqc, aes(
                x="sample_name", y=col, fill = "type")
                )
                + geom_col(stat='identity', poDuplicatessition='dodge')
                + labs(x="", y=col, title=col)
                + coord_flip()
                + theme_classic()
                + theme(
                    axis_text_x=element_text(rotation=-15, hjust=0.1)
                    )
            )
    output_name = "output/qc_fusions_barplots_{}.png".format(col)
    plot.save(output_name, height=10, width=10)