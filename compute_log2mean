
compute_log2mean <- function(data_list, replace_zero = FALSE, window_size = 10, cores = NA){

# Compute the vector of log2 mean ratios for each Hi-C row into list in parallel or not

        if (is.na(cores)){
                result <- lapply(X = data_list, FUN = function(data_table){
                        contact <- unname(unlist(data_table[,3]))
                        position <- unname(unlist(data_table[,2]))
                        log2mean <- unlist(lapply(
                                X = position, FUN = .estimate_ma_mean, 
                                data_set = contact, window_size = window_size))
                        if (replace_zero){log2mean[is.na(log2mean)] <- 0}
                        res <- cbind(data_table[,-3], log2mean)
                        res <- res[!is.nan(res[,3]),]
                        if (nrow(res) == 0){res <- rbind(res, c(NA, NA, NA))}
                        return(res)})
        } else {
                result <- parallel::mclapply(
                        1:length(data_list),
                        function(row_ind){
                                data_table <- data_list[[row_ind]]
                                contact <- unname(unlist(data_table[,3]))
                                position <- unname(unlist(data_table[,2]))
                                log2mean <- unlist(lapply(
                                        X = position, FUN = .estimate_ma_mean, 
                                        data_set = contact, window_size = window_size))
                                if (replace_zero){log2mean[is.na(log2mean)] <- 0}
                                res <- cbind(data_table[,-3], log2mean)
                                res <- res[!is.nan(res[,3]),]
                                if (nrow(res) == 0){res <- rbind(res, c(NA, NA, NA))}
                                return(res)},
                        mc.cores = cores)
                names(result) <- names(data_list)
        }
        
# Return the resulting matrix   

return(result)

}

# Compute log2 mean ratio at column position x based on the interactions within the given Hi-C row

.estimate_ma_mean <- function(x, data_set, window_size){
        if (!is.na(x)){
                if (x <= window_size){result <- NaN}
                if (x >= (length(data_set) - window_size + 1)){result <- NaN}
                if (x > window_size & x < (length(data_set) - window_size + 1)){
                        mean_left <- .estimate_single_ma_mean(
                                data_set = data_set, 
                                x_start = x - window_size, 
                                x_end = x - 1)
                        mean_right <- .estimate_single_ma_mean(
                                data_set = data_set, 
                                x_start = x + 1, 
                                x_end = x + window_size)
                        if (is.na(mean_left) | is.na(mean_right)){result <- NA}
                        if (!is.na(mean_left) & !is.na(mean_right)){result <- log2(mean_right/mean_left)}}
        } else {result <- NA}
        return(result)
}

# Compute MA mean estimate from provided interval based on the interactions within the given Hi-C row

.estimate_single_ma_mean <- function(data_set, x_start, x_end){
        select_set <- data_set[x_start:x_end]
        result <- mean(select_set, na.rm = TRUE)
        return(result)
}
