performOverlaps <- function(all.rearrs, bps.temp.gr, replis.gr, rfd.gr, genes.gr ) {
    # perform overlaps with SVs, replication timing, direction, and gene expression
    # all.rearrs: SVs in the bedpe format; row names are SV IDs
    # overlap with repliseq (average at 2 breakpoints)
    

    # overlaps between breakpoints and replication timing (genomic ranges objects)
    bps_repliseq.temp.df <- as.data.frame(findOverlaps(bps.temp.gr, replis.gr))
    
    bps_repliseq.temp.df$sv_id_global <- bps.temp.gr$sv_id_global[bps_repliseq.temp.df$queryHits]
    bps_repliseq.temp.df$repliseq <- replis.gr$V4[bps_repliseq.temp.df$subjectHits]

    bps_repliseq.temp.df_bp1 <- subset(bps_repliseq.temp.df, grepl('_bp1', bps_repliseq.temp.df$sv_id_global))
    bps_repliseq.temp.df_bp2 <- subset(bps_repliseq.temp.df, grepl('_bp2', bps_repliseq.temp.df$sv_id_global))
    bps_repliseq.temp.df_mp <- subset(bps_repliseq.temp.df, grepl('_midpoint', bps_repliseq.temp.df$sv_id_global))
    
    all.rearrs$repliseq_bp1 <- NA
    all.rearrs[gsub('_bp.*','',bps_repliseq.temp.df_bp1$sv_id_global),'repliseq_bp1'] <- bps_repliseq.temp.df_bp1$repliseq
    all.rearrs$repliseq_bp2 <- NA
    all.rearrs[gsub('_bp.*','',bps_repliseq.temp.df_bp2$sv_id_global),'repliseq_bp2'] <- bps_repliseq.temp.df_bp2$repliseq
    all.rearrs$repliseq_avg <- rowMeans(all.rearrs[,c('repliseq_bp1','repliseq_bp2')],na.rm=T)
    all.rearrs$repliseq_mp <- NA
    if (nrow(bps_repliseq.temp.df_mp)>0) {
        all.rearrs[gsub('_midpoint','',bps_repliseq.temp.df_mp$sv_id_global),'repliseq_mp'] <- bps_repliseq.temp.df_mp$repliseq
        all.rearrs$repliseq_mp[is.na(all.rearrs$repliseq_mp)] <- all.rearrs$repliseq_avg[is.na(all.rearrs$repliseq_mp)]
    }

    
    # overlap with replication fork direction
    rfd.gr$rf_dir <- NA
    rfd.gr$rf_dir[rfd.gr$score<0] <- 'L'
    rfd.gr$rf_dir[rfd.gr$score>0] <- 'R'
    
    bps_rfd.temp.df <- as.data.frame(findOverlaps(bps.temp.gr, rfd.gr))
    bps_rfd.temp.df$sv_id_global <- bps.temp.gr$sv_id_global[bps_rfd.temp.df$queryHits]
    bps_rfd.temp.df$rf_dir <- rfd.gr$rf_dir[bps_rfd.temp.df$subjectHits]
    bps_rfd.temp.df_bp1 <- subset(bps_rfd.temp.df, grepl('_bp1', bps_rfd.temp.df$sv_id_global))
    bps_rfd.temp.df_bp2 <- subset(bps_rfd.temp.df, grepl('_bp2', bps_rfd.temp.df$sv_id_global))
    
    bps_iz.temp.df <- as.data.frame(findOverlaps(bps.temp.gr, iz.gr))
    bps_iz.temp.df$sv_id_global <- bps.temp.gr$sv_id_global[bps_iz.temp.df$queryHits]
    bps_iz.temp.df_bp1 <- subset(bps_iz.temp.df, grepl('_bp1', bps_iz.temp.df$sv_id_global))
    bps_iz.temp.df_bp2 <- subset(bps_iz.temp.df, grepl('_bp2', bps_iz.temp.df$sv_id_global))
    
    
    all.rearrs$rf_dir_bp1 <- NA
    all.rearrs[gsub('_bp.*','',bps_rfd.temp.df_bp1$sv_id_global),'rf_dir_bp1'] <- bps_rfd.temp.df_bp1$rf_dir
    #all.rearrs[gsub('_bp.*','',bps_iz.temp.df_bp1$sv_id_global),'rf_dir_bp1'] <- 'IZ'
    all.rearrs$rf_dir_bp2 <- NA
    all.rearrs[gsub('_bp.*','',bps_rfd.temp.df_bp2$sv_id_global),'rf_dir_bp2'] <- bps_rfd.temp.df_bp2$rf_dir
    #all.rearrs[gsub('_bp.*','',bps_iz.temp.df_bp2$sv_id_global),'rf_dir_bp2'] <- 'IZ'
    all.rearrs$rf_dir <- paste(all.rearrs$rf_dir_bp1, all.rearrs$rf_dir_bp2, sep = "_")

    # overlap with genes
    # still need to account for overlap with multiple genes!
    bps_expr.temp.df <- as.data.frame(findOverlaps(bps.temp.gr, genes.gr))
    bps_expr.temp.df$sv_id_global <- bps.temp.gr$sv_id_global[bps_expr.temp.df$queryHits]
    bps_expr.temp.df$expr_fpkm <- genes.gr$expr_fpkm[bps_expr.temp.df$subjectHits]
    
    bps_expr.temp.df_bp1 <- subset(bps_expr.temp.df, grepl('_bp1', bps_expr.temp.df$sv_id_global))
    bps_expr.temp.df_bp2 <- subset(bps_expr.temp.df, grepl('_bp2', bps_expr.temp.df$sv_id_global))
    
    
    all.rearrs$expr_fpkm_bp1 <- 0
    all.rearrs[gsub('_bp.*','',bps_expr.temp.df_bp1$sv_id_global),'expr_fpkm_bp1'] <- bps_expr.temp.df_bp1$expr_fpkm
    all.rearrs$expr_fpkm_bp2 <- 0
    all.rearrs[gsub('_bp.*','',bps_expr.temp.df_bp2$sv_id_global),'expr_fpkm_bp2'] <- bps_expr.temp.df_bp2$expr_fpkm
    all.rearrs$expr_fpkm_avg <- rowMeans(all.rearrs[,c('expr_fpkm_bp1','expr_fpkm_bp2')],na.rm=T)

    result <- list()
    result[['all.rearrs']] <- all.rearrs

    # repliseq bins
    v <- c(all.rearrs$repliseq_bp1, all.rearrs$repliseq_bp2)
    quantiles <- quantile(v, probs = seq(0, 1, 0.2), na.rm = TRUE)
    labels <- c('very early', 'early', 'medium', 'late', 'very late')
    all.rearrs$repliseq_bp1_quantile <- cut(all.rearrs$repliseq_bp1, breaks = quantiles, include.lowest = TRUE, labels = labels)
    all.rearrs$repliseq_bp2_quantile <- cut(all.rearrs$repliseq_bp2, breaks = quantiles, include.lowest = TRUE, labels = labels)
    all.rearrs$repliseq_avg_quantile <- cut(all.rearrs$repliseq_avg, breaks = quantiles, include.lowest = TRUE, labels = labels)
    # this will be a 3D table
    t_repli <- table(all.rearrs$sample, all.rearrs$label, all.rearrs$repliseq_avg_quantile)
    result[['t_repli']] <- t_repli
    if (length(t_repli)>2) {
        result[['t_repli_2d']] <- collapseTable3D2D(t_repli)
    }
    # if any if the breakpoints have 
    if (nrow(bps_repliseq.temp.df_mp)>0) {
        all.rearrs$repliseq_mp_quantile <- cut(all.rearrs$repliseq_mp, breaks = quantiles, include.lowest = TRUE, labels = labels)
        t_repli_mp <- table(all.rearrs$sample, all.rearrs$label, all.rearrs$repliseq_mp_quantile)
        result[['t_repli_mp']] <- t_repli_mp
        if (length(t_repli_mp)>2) {
            result[['t_repli_mp_2d']] <- collapseTable3D2D(t_repli_mp)
        }
    }
    
    
    
    # expression bins
    v <- c(all.rearrs$expr_fpkm_bp1, all.rearrs$expr_fpkm_bp2)
    quantiles <- c(0,quantile(v[v>0], probs = seq(0, 1, 0.25), na.rm = TRUE))
    labels <- c('no_expr', 'low', 'medium low', 'medium high', 'high')
    all.rearrs$expr_bp1_quantile <- 'no_expr'
    all.rearrs$expr_bp2_quantile <- 'no_expr'
    all.rearrs$expr_avg_quantile <- 'no_expr'
    all.rearrs$expr_bp1_quantile <- cut(all.rearrs$expr_fpkm_bp1, breaks = quantiles, include.lowest = TRUE, labels = labels)
    all.rearrs$expr_bp2_quantile <- cut(all.rearrs$expr_fpkm_bp2, breaks = quantiles, include.lowest = TRUE, labels = labels)
    all.rearrs$expr_avg_quantile <- cut(all.rearrs$expr_fpkm_avg, breaks = quantiles, include.lowest = TRUE, labels = labels)
    # another 3D table
    t_expr <- table(all.rearrs$sample, all.rearrs$label, all.rearrs$expr_avg_quantile)
    result[['t_expr']] <- t_expr
    if (length(dim(t_expr))>2) {
        result[['t_expr_2d']] <- collapseTable3D2D(t_expr)
    }
        
    # RFD bins
    # another 3D table
    t_rfd <- table(all.rearrs$sample, all.rearrs$label, all.rearrs$rf_dir)
    t_rfd <- t_rfd[,,c('L_L', 'L_R', 'R_L', 'R_R')]
    result[['t_rfd']] <- t_rfd

    if (length(dim(t_rfd))>2) {
        result[['t_rfd_2d']] <- collapseTable3D2D(t_rfd)
    }
    
    return(result)
    
}

collapseTable3D2D <- function(table3D) {
    
    col_names_2 <- dimnames(table3D)[[2]] # Original second-dimension column names
    col_names_3 <- dimnames(table3D)[[3]]   # Original third-dimension names
    
    table2D <- matrix(table3D, nrow = dim(table3D)[1], ncol=length(col_names_2) * length(col_names_3), byrow=FALSE)
    new_col_names <- as.vector(outer(col_names_2, col_names_3, paste, sep = "_"))
    colnames(table2D) <- new_col_names
    rownames(table2D) <-  dimnames(table3D)[[1]]
    return(table2D)
}    

    