library(ggplot2)
message("Choose master table file")
master_table <- read.delim(file.choose())

add_rep_status <- function(master_table, rep_cutoff) {
    master_table$Repeat.Status <- factor(round(master_table$Repeat.Proportion
                                               - rep_cutoff + 0.5))
    return(master_table)
}

plot_base <- function(rep_cut) {
    if (!missing(rep_cut)) {
        master_table <- add_rep_status(master_table, rep_cut)
        p <- ggplot(data=master_table,
                aes(fill=Repeat.Status))
        p <- p +
             scale_fill_discrete(name=paste("Rep Status\n(cutoff=", rep_cut, ")", sep=""),
                                 breaks=c(0, 1),
                                 labels=c("Non-repeat", "Repeat"))
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
    matrix_scores <- master_table[!is.na(master_table$Max.Matrix.Score)
                                  & master_table$Max.Matrix.Score >= matrix_cut, ]
    return(nrow(matrix_scores))
}

get_num_rep_and_high_matrix <- function(rep_cut, matrix_cut) {
    matrix_scores <- master_table[!is.na(master_table$Max.Matrix.Score)
                                  & master_table$Max.Matrix.Score >= matrix_cut, ]
    reps <- master_table[master_table$Repeat.Proportion >= rep_cut, ]
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
         geom_point(aes(Max.MACS.Score, Max.Matrix.Score,
                        color=Repeat.Proportion)) +
         xlab("Max MACS Score") +
         ylab("Max Matrix Score") +
         scale_colour_gradient(name="% Repeats",
                               labels=paste(seq(0, 100, by=25),
                                            "%", sep=""))
    return(p)
}

plot_fe_vs_matrix <- function() {
    p <- plot_base()
    p <- p +
         geom_point(aes(Max.Fold.Enrichment, Max.Matrix.Score,
                        color=Repeat.Proportion)) +
         xlab("Max Fold Enrichment") +
         ylab("Max Matrix Score") +
         scale_colour_gradient(name="% Repeats",
                               labels=paste(seq(0, 100, by=25),
                                            "%", sep=""))
    return(p)
}

plot_macs_vs_fe <- function() {
    p <- plot_base()
    p <- p +
         geom_point(aes(Max.MACS.Score, Max.Fold.Enrichment,
                        color=Repeat.Proportion)) +
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

plot_chr_vs_macs <- function(rep_cut) {
    if (!missing(rep_cut)) {
        p <- plot_base(rep_cut)
    }
    else {
        p <- plot_base()
    }
    p <- p + geom_boxplot(aes(x=chr, y=Max.MACS.Score)) +
         xlab("Chromosome") +
         ylab("Max MACS Score")
    return(p)
}

plot_matrix_score_hist <- function(rep_cut) {
    if (!missing(rep_cut)) {
        p <- plot_base(rep_cut)
    }
    else {
        p <- plot_base()
    }
    p <- p +
         geom_histogram(aes(x=Max.Matrix.Score),
                        binwidth=.1) +
         xlab("Max Matrix Score")
    return(p)
}

plot_anno_cat_hist <- function(rep_cut) {
    if (!missing(rep_cut)) {
        p <- plot_base(rep_cut)
    }
    else {
        p <- plot_base()
    }
    p <- p +
         geom_bar(aes(x=reorder(Annotation.Category,
                                Annotation.Category,
                                function(x)-length(x)))) +
         xlab("Annotation Category") +
         theme(axis.text.x=element_text(angle=45, hjust=1))
    return(p)
}

plot_anno_cat_vs_rep <- function() {
    p <- plot_base() +
         geom_point(aes(x=reorder(Annotation.Category,
                                  Annotation.Category,
                                  function(x)-length(x)),
                        y=Repeat.Proportion)) +
         xlab("Annotation Category") +
         ylab("% Repeats") +
         theme(axis.text.x=element_text(angle=45, hjust=1)) +
         scale_y_continuous(labels=scales::percent)
    return(p)
}

plot_num_repeat_vs_length <- function() {
    p <- plot_base() +
         geom_point(aes(x=Repeat.Count, y=Length,
                        color=Repeat.Proportion))
    return(p)
}
