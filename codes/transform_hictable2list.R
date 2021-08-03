
transform_hictable2list <- function(data_table, replace_zero = FALSE, limit_size = NA, window_size = NA, cores = NA){

# Find gaps between consevutive loci

        locus_start <- sort(unique(c(
                paste0(data_table[,1], "_", data_table[,2]),
                paste0(data_table[,4], "_", data_table[,5]))))
        locus_end <- sort(unique(c(
                paste0(data_table[,1], "_", data_table[,3]),
                paste0(data_table[,4], "_", data_table[,6]))))
        gaps <- cbind(locus_start[-1], locus_end[-length(locus_end)])
        gaps <- gaps[
                gaps[,1] != gaps[,2] &
                as.character(unlist(strsplit(gaps[,1], "_"))[c(T,F)]) == as.character(unlist(strsplit(gaps[,1], "_"))[c(T,F)]),]
        bins_full <- rbind(cbind(locus_start, locus_end), gaps)
        locus <- paste0(bins_full[,1], "_", as.character(unlist(strsplit(bins_full[,2], "_"))[c(F,T)]))
        bins_full <- cbind(bins_full, locus)
        bins_full <- bins_full[order(bins_full[,1]),]
        bins_n <- nrow(bins_full)
        
# Create an index for each pair of interacting loci

        contact_table <- data.frame(
                "index" = paste0(
                        match(paste0(data_table[,4], "_", data_table[,5], "_", data_table[,6]), bins_full[,3]), "_",
                        match(paste0(data_table[,1], "_", data_table[,2], "_", data_table[,3]), bins_full[,3])),
                "contact" = data_table[,7])
        contact_table[,1] <- as.character(contact_table[,1])         
        
# Specify column positions to extract for each Hi-C row depending on full/restricted triangle  
        
        if (is.na(window_size)){window_size <- 0}
        if (is.na(limit_size)){
                col_start <- c(
                        rep(NA, times = 2*window_size + 1),
                        rep(1, times = bins_n - 2*window_size - 1))
                col_end <- c(
                        rep(NA, times = 2*window_size + 1),
                        (2*window_size + 1):(bins_n - 1))
        } else {
                col_start <- c(
                        rep(NA, times = 2*window_size + 1),
                        rep(1, times = limit_size - 1),
                        1:(bins_n - 2*window_size - limit_size))
                col_end <- c(
                        rep(NA, times = 2*window_size + 1),
                        (2*window_size + 1):(bins_n - limit_size),
                        rep(bins_n - limit_size, times = limit_size - 1))}
                        
# Extract interactions column-wise for each Hi-C row into list in parallel or not including lost interactions (as NAs or zeros)

        result <- split(bins_full[,3], bins_full[,3])
        if (is.na(cores)){
                result <- lapply(X = result, FUN = function(row_locus){
                        row_ind <- match(row_locus, bins_full[,3])
                        if (is.na(col_start[row_ind])){
                                res <- data.frame("col_locus" = NA, "col_index" = NA, "contact" = NA)
                        } else {
                                col_ind <- col_start[row_ind]:col_end[row_ind]
                                loci_ind <- paste0(row_ind, "_", col_ind)
                                res <- data.frame(
                                        "col_locus" = bins_full[col_ind, 3],
                                        "col_ind" = col_ind,
                                        "contact" = contact_table[match(loci_ind, contact_table[,1]),2])
                                if (replace_zero){if (length(is.na(res[,3])) != 0){res[is.na(res[,3]), 3] <- 0}}}
                        res[,1] <- as.character(res[,1])
                        return(res)})
        } else {
                result <- parallel::mclapply(
                        1:length(result),
                        function(row_ind){
                                if (is.na(col_start[row_ind])){
                                        res <- data.frame("col_locus" = NA, "col_index" = NA, "contact" = NA)
                                } else {
                                        col_ind <- col_start[row_ind]:col_end[row_ind]
                                        loci_ind <- paste0(row_ind, "_", col_ind)
                                        res <- data.frame(
                                                "col_locus" = bins_full[col_ind, 3],
                                                "col_ind" = col_ind,
                                                "contact" = contact_table[match(loci_ind, contact_table[,1]),2])
                                        if (replace_zero){if (length(is.na(res[,3])) != 0){res[is.na(res[,3]), 3] <- 0}}}
                                res[,1] <- as.character(res[,1])
                                return(res)},
                        mc.cores = cores)
                names(result) <- bins_full[,3]}
                
# Return the resulting list

        return(result)
        
}
