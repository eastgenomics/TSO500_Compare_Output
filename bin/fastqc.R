setwd("/home/aisha/Documents/Projects/cancer/pancan/validation_run2/compare_data/rna_fastqc_to_compare/")
library(reshape2)
library(ggplot2)

dat <- read.csv("fastqc_validations_all.txt", sep = "\t")
char <- strsplit(dat$Sample.Name,"-")
dat$sample_names <- sapply(char,"[[",1)
dat$batch <- sapply(char,"[",3)
dat$type <- rep(NA,nrow(dat))

for (i in 1:nrow(dat)){
  if (grepl("_", dat[i,"sample_names"], fixed = TRUE)){
    dat[i,"type"] <- "Neat 22"
  }
  else if (endsWith(dat[i,"sample_names"], 'R')){
    dat[i,"type"] <- "1:4 dilution 23PCR1"
  }
  else {
    dat[i,"type"] <- "1:2 dilution 23PCR1"
  }
}
# anything 24PCR1 is 1:2 dilution 24PCR1
for (i in 1:nrow(dat)){
  if (grepl("24PCR1", dat[i,"batch"], fixed = TRUE)){
    dat[i,"type"] <- "1:2 dilution 24PCR1"
  }
}
# anything that was done on TSO is 1:2 dilution 24PCR1
for (i in 1:nrow(dat)){
  if (grepl("23TSO", dat[i,"batch"], fixed = TRUE)){
    dat[i,"type"] <- "TSO500 RNA"
  }
}

dat$sample_names <- sapply(strsplit(dat$sample_names,"_"),"[[",1)
dat$sample_names <- sapply(strsplit(dat$sample_names,"R"),"[[",1)
# remove the lane info from samplename 
dat$sample_names <- gsub("\\_.*","", sapply(char,"[[",1))
# remove the "R" for repeat
char <- strsplit(dat$sample_names,"-")
dat$sample_names <- gsub("\\R.*","", sapply(char,"[[",1))

plot_metric <- function(dat, metric) { 
  dat_long <- melt(dat,
                   id.vars=c("sample_names", "type"),
                   measure.vars=c(metric),
                   variable.name="Metric",
                   value.name=metric
  )
  dat_long$type <- factor(dat_long$type, levels = c("Neat 22","1:2 dilution 23PCR1", "1:2 dilution 24PCR1", "1:4 dilution 23PCR1", "TSO500 RNA"))
  
  ggplot(dat_long, aes(x=sample_names, y=get(metric), fill=type)) +
    geom_bar(stat='identity', position='dodge') +
    theme_classic() +
    theme(axis.text.x=element_text(size=15, angle=90,hjust=0.95,vjust=0.2)) +
    labs(title = metric, x = "", y = metric) +
    theme(plot.title = element_text(hjust = 0.5))
  
}
metrics = c("Dups", "GC", "M.Seqs")
pdf("/home/aisha/Documents/Projects/cancer/pancan/validation_run2/compare_data/output/fastqc_validation2_plots.pdf",
    width = 15,
    height = 10)
for (metric in metrics){
  print(plot_metric(dat,metric))
}
dev.off()

