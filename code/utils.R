
#' Add chromosome number and SNP positions for SNPs started by Chr*
#' 
#' @param snpdf Data frame of SNP positions with the mandatory columns
#' 'SNP' (SNP name), 'chr' (chromosome number), and 'pos' (SNP position).
#' @return The same data frame as in the input, but updated with the columns
#' 'chr', 'start' and 'end'  
#' @noRd
add_coordinates <- function(snpdf = NULL) {
    # Add chromosome info
    snpdf <- snpdf[!startsWith(snpdf$SNP, "BARC"), ]
    snpdf$chr[is.na(snpdf$chr)] <- sapply(
        strsplit(snpdf$SNP[is.na(snpdf$chr)], "-"), `[`, 1
    )
    
    # Add start and end positions
    snpdf$start <- snpdf$pos
    snpdf$start[is.na(snpdf$start)] <- sapply(
        strsplit(snpdf$SNP[is.na(snpdf$start)], "-"),
        `[`, 2
    ) 
    
    snpdf$end <- snpdf$pos
    snpdf$end[is.na(snpdf$end)] <- sapply(
        strsplit(snpdf$SNP[is.na(snpdf$end)], "-"), 
        tail, 1
    )
    return(snpdf)
}


#' Left-pad chromosome 1-9 with zeroes if necessary
#' 
#' @param snpdf Data frame of SNPs and positions. Column 'chr' is mandatory.
#' 
#' @return The same input data frame, but with fixed names.
#' @noRd
chr_rename <- function(snpdf) {
    snpdf$chr <- gsub("Chr1$", "Chr01", snpdf$chr)
    snpdf$chr <- gsub("Chr2$", "Chr02", snpdf$chr)
    snpdf$chr <- gsub("Chr3$", "Chr03", snpdf$chr)
    snpdf$chr <- gsub("Chr4$", "Chr04", snpdf$chr)
    snpdf$chr <- gsub("Chr5$", "Chr05", snpdf$chr)
    snpdf$chr <- gsub("Chr6$", "Chr06", snpdf$chr)
    snpdf$chr <- gsub("Chr7$", "Chr07", snpdf$chr)
    snpdf$chr <- gsub("Chr8$", "Chr08", snpdf$chr)
    snpdf$chr <- gsub("Chr9$", "Chr09", snpdf$chr)
    
    return(snpdf)
}


#' Get workspace size
#' 
#' @return Messages.
#' @noRd
workspace_size <- function(env = globalenv()) {
    env_vars <- ls(envir = env)
    size <- 0
    for (x in env_vars) {
        thisSize <- object.size(get(x))
        size <- size + thisSize
        message(x, " = ", appendLF = F); print(thisSize, units='auto')
    }
    message("Total workspace is ",appendLF = F); print(size, units='auto')
}