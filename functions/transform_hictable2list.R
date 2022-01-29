transform_hictable2list <- function(data_table, direction, fill_empty = FALSE, replace_zero = FALSE, limit_size = NA, resolution = NA, cores = NA){
# Remove inter-chromosomal interactions
        data_table <- data_table[data_table[,1] == data_table[,4],]
        colnames(data_table) <- c("chr", "position", "position", "chr", "position", "position", "contact")
# Find gaps between consevutive loci
        if (fill_empty){
                bins_full <- rbind(data_table[,c(1,2)], data_table[,c(1,3)], data_table[,c(4,5)], data_table[,c(4,6)])
                unique_ind <- which(!duplicated(paste0(bins_full[,1], "_", bins_full[,2])))
                bins_full <- bins_full[unique_ind,]
                order_full <- order(as.character(bins_full[,1]), as.numeric(bins_full[,2]))
                bins_full <- bins_full[order_full,]
                bins_full[,1] <- as.character(bins_full[,1])
                chr_end <- lapply(X = split(bins_full[,2], bins_full[,1]), FUN = max)
                chr_name <- as.character(names(chr_end))
                chr_end <- as.numeric(unname(unlist(chr_end)))
                if (is.na(resolution)){
                        for (i in 1:length(chr_name)){
                                remove_ind <- which(bins_full[,1] == chr_name[i] & bins_full[,2] == chr_end[i])
                                bins_full <- bins_full[-remove_ind,]}
                } else {
                        bins_full <- c(NA,NA)
                        for (i in 1:length(chr_name)){
                                if (chr_end[i] %% resolution == 0){
                                        bins_full <- rbind(bins_full, cbind(
                                                rep(chr_name[i], times = chr_end[i]/resolution),
                                                seq(from = 0, to = chr_end[i] - resolution, by = resolution)))
                                } else {
                                        bins_full <- as.data.frame(rbind(bins_full, cbind(
                                                rep(chr_name[i], times = chr_end[i]%/%resolution + 1),
                                                seq(from = 0, to = chr_end[i], by = resolution))))}}
                            bins_full <- bins_full[-1,]
                            bins_full[,1] <- as.character(bins_full[,1])
                            bins_full[,2] <- as.numeric(as.character(bins_full[,2]))}
        }else{
                bins_full <- rbind(data_table[,c(1,2)], data_table[,c(4,5)])
                unique_ind <- which(!duplicated(paste0(bins_full[,1], "_", bins_full[,2])))
                bins_full <- bins_full[unique_ind,]
                order_full <- order(as.character(bins_full[,1]), as.numeric(bins_full[,2]))
                bins_full <- bins_full[order_full,]
                bins_full[,1] <- as.character(bins_full[,1])}
        bins_n <- nrow(bins_full)
        bins_locus <- paste0(bins_full[,1], "_", bins_full[,2])
# Create an index for each pair of interacting loci (lower triangle)
        contact_table <- data.frame(
                "row_ind" = match(paste0(data_table[,4], "_", data_table[,5]), bins_locus),
                "col_ind" = match(paste0(data_table[,1], "_", data_table[,2]), bins_locus),
                "contact" = data_table[,7])
        contact_table$coordinate <- as.character(paste0(contact_table$row_ind, "-", contact_table$col_ind))
        contact_table <- contact_table[,c(4,3)]
# Specify indeces to extract for each Hi-C row/column depending on full/restricted triangle (Fig. 6A)   
        if (is.na(limit_size)){
                if (direction == "row"){
                        ind_start <- c(NA, rep(1, times = bins_n - 1))
                        ind_end <- c(NA, 1:(bins_n - 1))}
                if (direction == "column"){
                        ind_start <- c(2:bins_n, NA)
                        ind_end <- c(rep(bins_n, times = bins_n - 1), NA)}
        } else {
                if (limit_size + 1 >= bins_n){
                        warning(paste0("The parameter limit_size exceeds the number of bins available. Non-restricted matrix is used."))
                        if (direction == "row"){
                                ind_start <- c(NA, rep(1, times = bins_n - 1))
                                ind_end <- c(NA, 1:(bins_n - 1))}
                        if (direction == "column"){
                                ind_start <- c(2:bins_n, NA)
                                ind_end <- c(rep(bins_n, times = bins_n - 1), NA)}
                } else {
                        if (direction == "row"){
                                ind_start <- c(NA, rep(1, times = limit_size - 1), 1:(bins_n - limit_size))
                                ind_end <- c(NA, 1:(bins_n - limit_size), rep(bins_n - limit_size, times = limit_size - 1))}
                        if (direction == "column"){
                                ind_start <- c(rep(limit_size + 1, times = limit_size - 1), (limit_size + 1):(bins_n), NA)
                                ind_end <- c((limit_size + 1):bins_n, rep(bins_n, times = limit_size - 1), NA)}}}
# Extract interactions within each Hi-C row/column into list in parallel or not including lost interactions (as NAs or zeros)
        if (is.na(cores)){
                result <- lapply(X = 1:length(bins_locus), FUN = function(locus_ind){
                        if (is.na(ind_start[locus_ind])){
                                res <- data.frame("V1" = NA, "V2" = NA, "V3" = NA)
                        } else {
                                position_ind <- ind_start[locus_ind]:ind_end[locus_ind]
                                if (direction == "row"){coords <- as.character(paste0(locus_ind, "-", position_ind))}
                                if (direction == "column"){coords <- as.character(paste0(position_ind, "-", locus_ind))}
                                res <- data.frame(
                                        "V1" = bins_locus[position_ind],
                                        "V2" = position_ind,
                                        "V3" = contact_table[match(coords, contact_table[,1]),2])
                                if (replace_zero){if (length(is.na(res[,3])) != 0){res[is.na(res[,3]), 3] <- 0}}}
                        if (direction == "row"){colnames(res) <- c("col_locus", "col_index", "contact")}
                        if (direction == "column"){colnames(res) <- c("row_locus", "row_index", "contact")}
                        res[,1] <- as.character(res[,1])
                        return(res)})
        } else {
                result <- parallel::mclapply(
                        1:length(bins_locus),
                        function(locus_ind){
                                if (is.na(ind_start[locus_ind])){
                                        res <- data.frame("V1" = NA, "V2" = NA, "V3" = NA)
                                } else {
                                        position_ind <- ind_start[locus_ind]:ind_end[locus_ind]
                                        if (direction == "row"){coords <- paste0(locus_ind, "-", position_ind)}
                                        if (direction == "column"){coords <- paste0(position_ind, "-", locus_ind)}
                                        res <- data.frame(
                                                "V1" = bins_locus[position_ind],
                                                "V2" = position_ind,
                                                "V3" = contact_table[match(coords, contact_table[,1]),2])
                                        if (replace_zero){if (length(is.na(res[,3])) != 0){res[is.na(res[,3]), 3] <- 0}}}
                                if (direction == "row"){colnames(res) <- c("col_locus", "col_index", "contact")}
                                if (direction == "column"){colnames(res) <- c("row_locus", "row_index", "contact")}
                                res[,1] <- as.character(res[,1])
                                return(res)},
                        mc.cores = cores)
                names(result) <- bins_locus}
# Return the resulting list
        return(result)
}
