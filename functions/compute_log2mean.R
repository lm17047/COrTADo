compute_log2mean <- function(data_list, replace_zero = FALSE, window_size, cores = NA){
# Compute the vector of log2 mean ratios for each Hi-C row/column into list in parallel or not
        if (is.na(cores)){
                result <- lapply(X = data_list, FUN = function(data_table){
                        contact <- unname(unlist(data_table[,3]))
                        position <- unname(unlist(data_table[,2]))
                        log2mean <- unlist(lapply(
                                X = position, FUN = .estimate_ma_mean, 
                                data_set = contact, data_index = position,
                                window_size = window_size,
                                start_limit = position[1] + window_size - 1, 
                                end_limit = position[length(position)] - window_size + 1)) 
                        if (replace_zero){log2mean[is.na(log2mean)] <- 0}
                        res <- cbind(data_table[,-3], log2mean)
                        res <- res[!is.nan(res[,3]),]
                        res_names <- colnames(res)
                        if (nrow(res) == 0){
                                res <- rbind(res, c(NA, NA, NA))
                                colnames(res) <- res_names}
                        return(res)})
        } else {
                result <- parallel::mclapply(
                        1:length(data_list),
                        function(list_ind){
                                data_table <- data_list[[list_ind]]
                                contact <- unname(unlist(data_table[,3]))
                                position <- unname(unlist(data_table[,2]))
                                log2mean <- unlist(lapply(
                                        X = position, FUN = .estimate_ma_mean, 
                                        data_set = contact, data_index = position,
                                        window_size = window_size,
                                        start_limit = position[1] + window_size - 1,
                                        end_limit = position[length(position)] - window_size + 1)) 
                                if (replace_zero){log2mean[is.na(log2mean)] <- 0}
                                res <- cbind(data_table[,-3], log2mean)
                                res <- res[!is.nan(res[,3]),]
                                res_names <- colnames(res)
                                if (nrow(res) == 0){
                                        res <- rbind(res, c(NA, NA, NA))
                                        colnames(res) <- res_names}
                                return(res)},
                        mc.cores = cores)
                names(result) <- names(data_list)}
# Return the resulting matrix                                        
        return(result)
}
