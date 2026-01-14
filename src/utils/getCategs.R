getCategs <- function(order4) {
    order1 <- c('clustered', 'non-clustered')
    order2 <- c('del', 'tds', 'inv', 'trans')
    order3 <- c('1-10Kb', '10-100Kb', '100Kb-1Mb', '1Mb-10Mb','>10Mb')
    totalCats <- length(order1) * length(order2) * length(order3) * length(order4)   
    categs <- paste0( rep(order1, each=totalCats/length(order1)), '_', 
                     rep(rep(order2, each=length(order3)*length(order4)), length(order1)), '_',
                    rep(rep(order3, each=length(order4)), length(order1) * length(order2)), '_',
                    rep(order4, length(order1) * length(order2) * length(order3))
                    )
    categs[grepl( 'trans_.*b_',categs)] <- gsub('trans_.*b_','trans_',categs[grepl( 'trans_.*b_',categs)])
    categs <- categs[!duplicated(categs)]
    return(categs)
}