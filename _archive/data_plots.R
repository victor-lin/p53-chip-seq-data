library(ggplot2)
library(scales)
message("Choose master table file")
master_table <- read.delim(file.choose())


get_rep_interval_names <- function(rep_cutoffs) {
    # ex.
    # rep_cutoffs=c(0, .1, .9, 1)
    # return c("0-10%", "10%-90%", "90%-100%")
    rep_cutoffs <- percent(rep_cutoffs)
    return (sapply(1:(length(rep_cutoffs) - 1),
                   function(x) paste(rep_cutoffs[x], '-',
                                     rep_cutoffs[x + 1], sep='')))
}

add_rep_status <- function(master_table, rep_cutoffs) {
    rep_cutoffs <- c(0, rep_cutoffs, 1)
    interval_names <- get_rep_interval_names(rep_cutoffs)
    repeat_status <- interval_names[findInterval(master_table$repeat_proportion,
                                                 rep_cutoffs,
                                                 rightmost.closed=TRUE)]
    master_table$repeat_status <- factor(repeat_status)
    return(master_table)
}

plot_base <- function(rep_cutoffs) {
    if (!missing(rep_cutoffs)) {
        master_table <- add_rep_status(master_table, rep_cutoffs)
        p <- ggplot(data=master_table,
                aes(fill=repeat_status))
        p <- p +
             scale_fill_discrete(name="Repeat Status")
    }
    else {
        p <- ggplot(data=master_table)
    }
    p <- p +
         theme_bw() +
         theme(legend.key = element_blank())
    return(p)
}

get_num_high_matrix <- function(matrix_cut) {
    matrix_scores <- master_table[!is.na(master_table$P53match_score_max)
                                  & master_table$P53match_score_max >= matrix_cut, ]
    return(nrow(matrix_scores))
}

get_num_rep_and_high_matrix <- function(rep_cutoffs, matrix_cut) {
    matrix_scores <- master_table[!is.na(master_table$P53match_score_max)
                                  & master_table$P53match_score_max >= matrix_cut, ]
    reps <- master_table[master_table$repeat_proportion >= rep_cutoffs, ]
    rep_and_matrix <- merge(matrix_scores, reps)
    return(nrow(rep_and_matrix))
}

plot_rep_cutoffs <- function(matrix_cut, interval) {
    rep_seq <- seq(0, 1, by=interval)
    df <- data.frame(rep_cutoff=rep_seq)
    df['num_high_matrix'] = get_num_high_matrix(matrix_cut)
    df['num_rep'] <- sapply(df$rep_cutoff,
                            FUN=get_num_rep_and_high_matrix,
                            matrix_cut=matrix_cut)
    df['num_nonrep'] <- df$num_high_matrix - df$num_rep
    # reshape data
    library(reshape2)
    meltd <- melt(df, id.vars=1:2)
    p <- ggplot(data=meltd, aes(x=rep_cutoff, y=value, fill=variable)) +
         geom_bar(stat="identity") +
         scale_x_continuous(breaks=rep_seq) +
         scale_y_continuous(breaks=seq(0, max(df$num_high_matrix), by=6)) +
         scale_fill_discrete(name="Rep Status",
                             breaks=c("num_rep", "num_nonrep"),
                             labels=c("Rep", "Non-rep")) +
         xlab("Rep cutoff") +
         ylab("# high matrix score") +
         scale_x_continuous(breaks=rep_seq, labels=scales::percent)
         theme_bw() +
         theme(legend.key=element_blank())
    return(p)
}

plot_macs_vs_matrix <- function() {
    p <- plot_base()
    p <- p +
         geom_point(aes(MACS_max, P53match_score_max,
                        color=repeat_proportion)) +
         xlab("Max MACS Score") +
         ylab("Max P53 Match Score") +
         scale_colour_gradient(name="% Repeats",
                               labels=paste(seq(0, 100, by=25),
                                            "%", sep=""))
    return(p)
}

plot_fe_vs_matrix <- function() {
    p <- plot_base()
    p <- p +
         geom_point(aes(FE_max, P53match_score_max,
                        color=repeat_proportion)) +
         xlab("Max Fold Enrichment") +
         ylab("Max P53 Match Score") +
         scale_colour_gradient(name="% Repeats",
                               labels=paste(seq(0, 100, by=25),
                                            "%", sep=""))
    return(p)
}

plot_macs_vs_fe <- function() {
    p <- plot_base()
    p <- p +
         geom_point(aes(MACS_max, FE_max,
                        color=repeat_proportion)) +
         xlab("Max MACS Score") +
         ylab("Max Fold Enrichment") +
         scale_colour_gradient(name="% Repeats",
                               labels=paste(seq(0, 100, by=25),
                                            "%", sep=""))
    return(p)
}

# plot_sample_peak_lengths <- function(rep) {
#     p <- ggplot(data=valid, aes(x=Sample_Name, y=sample_length)) +
#          geom_boxplot() +
#          theme_bw() +
#          theme(axis.text.x=element_text(angle=45, hjust=1)) +
#          xlab("Sample Name") +
#          ylab("Peak Length")
#     if (rep) {
#         p = p + geom_boxplot(aes(fill=Rep))
#     }
#     return(p)
# }

plot_chr_vs_macs <- function(rep_cutoffs) {
    if (!missing(rep_cutoffs)) {
        p <- plot_base(rep_cutoffs)
    }
    else {
        p <- plot_base()
    }
    p <- p + geom_boxplot(aes(x=chr, y=MACS_max)) +
         xlab("Chromosome") +
         ylab("Max MACS Score")
    return(p)
}

plot_matrix_score_hist <- function(rep_cutoffs) {
    if (!missing(rep_cutoffs)) {
        p <- plot_base(rep_cutoffs)
    }
    else {
        p <- plot_base()
    }
    p <- p +
         geom_histogram(aes(x=P53match_score_max),
                        binwidth=.1) +
         xlab("Max P53 Match Score")
    return(p)
}

plot_macs_score_hist <- function(rep_cutoffs) {
    if (!missing(rep_cutoffs)) {
        p <- plot_base(rep_cutoffs)
    }
    else {
        p <- plot_base()
    }
    p <- p +
         geom_histogram(aes(x=MACS_max), binwidth=100) +
         xlab("Max MACS Score") +
         scale_x_continuous(limits=c(50, NA))
    return(p)
}

plot_rep_count_hist <- function() {
    p <- plot_base() +
         geom_histogram(aes(x=repeat_count), binwidth=50) +
         xlab("Repeat Count") +
         ylab("Frequency")
    return(p)
}

plot_rep_percent_hist <- function(breaks) {
    # args:
    #   breaks: intervals for repeat percent.
    #           ex. c(0,1,50,90,100)
    p <- plot_base()
    if (missing(breaks)) {
        p <- p + geom_histogram(aes(x=repeat_proportion), binwidth=.01)
    }
    else {
        p <- p + geom_histogram(aes(x=repeat_proportion),
                                breaks=breaks / 100, fill="white", color="black")
    }
    p <- p +
         scale_x_continuous(labels=scales::percent) +
         xlab("Repeat Percent") +
         ylab("Frequency")
    return(p)
}

plot_anno_cat_hist <- function(rep_cutoffs) {
    if (!missing(rep_cutoffs)) {
        p <- plot_base(rep_cutoffs)
    }
    else {
        p <- plot_base()
    }
    p <- p +
         geom_bar(aes(x=reorder(annotation, annotation,
                                function(x)-length(x)))) +
         xlab("Annotation Category") +
         theme(axis.text.x=element_text(angle=45, hjust=1))
    return(p)
}

plot_anno_cat_vs_rep <- function() {
    p <- plot_base() +
         geom_point(aes(x=reorder(annotation, annotation,
                                  function(x)-length(x)),
                        y=repeat_proportion)) +
         xlab("Annotation Category") +
         ylab("% Repeats") +
         theme(axis.text.x=element_text(angle=45, hjust=1)) +
         scale_y_continuous(labels=scales::percent)
    return(p)
}

plot_rep_count_vs_length <- function() {
    p <- plot_base() +
         geom_point(aes(x=repeat_count, y=peak_length,
                        color=repeat_proportion)) +
         xlab("Repeat Count") +
         ylab("Peak Length")
    return(p)
}

plot_rep_count_vs_length_box <- function() {
    p <- plot_base() +
         geom_boxplot(aes(x=repeat_count, y=peak_length,
                          group=cut_width(repeat_count,
                                          width=500,
                                          boundary=1,
                                          closed="left"))) +
         xlab("Repeat Count") +
         ylab("Peak Length")
    return(p)
}

plot_sample_anno_distribution <- function() {
    # Stacked barplot of annotation distributions for samples.
    # Annotation stacks sorted by frequency in master table.
    message("Choose summary annotation file")
    anno_summary <- read.delim(file.choose())
    require(RColorBrewer)
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
    category.order <- sort(table(master_table$annotation),
                           decreasing=TRUE)
    anno_summary$annotation <- factor(anno_summary$annotation,
                                      levels=names(category.order))
    p <- ggplot(data=anno_summary) +
         geom_bar(aes(x=Sample.Name, fill=annotation),
                  position = "fill") +
         theme(axis.text.x=element_text(angle=45, hjust=1)) +
         scale_fill_manual(values=getPalette(length(category.order)),
                           name="Annotation Category") +
         scale_y_continuous(labels=scales::percent) +
         xlab("Sample Name") +
         ylab("Frequency")
    return(p)
}
