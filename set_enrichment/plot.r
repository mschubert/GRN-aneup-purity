library(cowplot)

odds = function(df)
    ggplot(df, aes(x=edges, y=estimate, color=method)) +
        geom_hline(yintercept=1, color="black", linetype="dashed") +
        geom_line() +
        geom_point() +
        facet_grid(expr ~ cohort) +
        scale_x_log10() +
        scale_y_log10(breaks=c(0.1,1,3,10,30,100)) +
        theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
              axis.text.x = element_text(angle=45, hjust=1)) +
        ggtitle("Odds ratio for links in the same CNA")

fp_total = function(df)
    ggplot(df, aes(x=edges, y=pmax(0,exp_fp), color=method)) +
        geom_line() +
        geom_point() +
        facet_grid(expr ~ cohort) +
        scale_x_log10() +
        scale_y_log10() +
        theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
              axis.text.x = element_text(angle=45, hjust=1)) +
        ggtitle("Expected number of false positive links")

fpr_genome = function(df)
    ggplot(df, aes(x=edges, y=fpr_total*100, color=method)) +
        geom_line() +
        geom_point() +
        facet_grid(expr ~ cohort, scales="free_y") +
        scale_x_log10() +
        theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
              axis.text.x = element_text(angle=45, hjust=1)) +
        ggtitle("FPR genome-wide")

fpr_segment = function(df)
    ggplot(df, aes(x=edges, y=pmax(0,fpr_seg*100), color=method)) +
        geom_line() +
        geom_point() +
        facet_grid(expr ~ cohort) +
        scale_x_log10() +
        theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
              axis.text.x = element_text(angle=45, hjust=1)) +
        ggtitle("FPR within CNA")
