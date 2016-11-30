library("optparse")
option_list = list(
    make_option("--bed_directory", type="character", default=NULL,
                help="BED files directory", metavar="character"),
    make_option(c("-o", "--out_file"), type="character", default=NULL,
                help="output filepath", metavar="character"),
    make_option("--fe_directory", type="character", default=NULL,
                help="FE files directory", metavar="character"),
    make_option("--ignore", type="character", default=NULL,
                help="ignore a chromosome from analysis", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

# TODO: check valid directory paths

# summary sample BED dataframe
setwd(args$bed_directory)
for (fname in list.files()) {
    sample.name <- unlist(strsplit(fname, "_peaks.bed"))
    sample.df <- read.delim(fname, header=FALSE,
                            col.names=c("chr", "start", "end", "PeakID",
                                        "MACS"))
    sample.df$Sample_Name <- sample.name
    if (fname == list.files()[1]) {
        summary.df <- sample.df
    }
    else {
        summary.df <- rbind(summary.df, sample.df)
    }
}

# FE dataframe
setwd(args$fe_directory)
for (fname in list.files()) {
    sample.name <- unlist(strsplit(fname, "_peaks.xls"))
    sample.df <- read.delim(fname, header=TRUE, skip=23)
    colnames(sample.df)[colnames(sample.df) == "fold_enrichment"] <- "FE"
    sample.df$Sample_Name <- sample.name
    if (fname == list.files()[1]) {
        fe.df <- sample.df
    }
    else {
        fe.df <- rbind(fe.df, sample.df)
    }
}

# some annotation script had offset `start` by 1, undo here
fe.df$start <- fe.df$start - 1

df <- merge(summary.df[, !(names(summary.df) %in% "PeakID")],
            fe.df[, c("chr", "start", "end",
                      "Sample_Name", "FE")],
            by=(c("chr", "start", "end", "Sample_Name")))

# add length column, reorder columns
df$length <- df$end - df$start
df <- df[, c("chr", "start", "end", "length",
             "Sample_Name", "MACS", "FE")]

if (!is.null(args$ignore)) {
    df <- df[df$chr != args$ignore, ]
}
write.table(df, file=args$out_file, sep="\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)
