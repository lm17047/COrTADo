window_size_diagnostics <- function(data_list, limit_size = NA, cores = NA){
        if (is.na(cores)){
                consecutive_na <- lapply(X = data_list, FUN = function(data_table){
                        if (nrow(data_table) > 1){
                                condition_na <- rep("F", times = nrow(data_table))
                                condition_na[is.na(data_table[,3])] <- "T"
                                condition_na <- paste0(condition_na[-length(condition_na)], condition_na[-1])
                                start_na <- which(condition_na == "FT") + 1
                                if (condition_na[1] == "TT" | 
                                    condition_na[1] == "TF"){
                                            start_na <- c(1, start_na)}
                                end_na <- which(condition_na == "TF")
                                if (condition_na[length(condition_na)] == "TT" | 
                                    condition_na[length(condition_na)] == "FT"){
                                                end_na <- c(end_na, length(condition_na) + 1)}
                                if (length(start_na) > 0){res <- end_na - start_na + 1} else {res <- NA}
                        } else {res <- NA}
                        return(res)})
        } else {
                consecutive_na <- parallel::mclapply(
                        1:length(data_list),
                        function(list_ind){
                                data_table <- data_list[[list_ind]]
                                if (nrow(data_table) > 1){
                                        condition_na <- rep("F", times = nrow(data_table))
                                        condition_na[is.na(data_table[,3])] <- "T"
                                        condition_na <- paste0(condition_na[-length(condition_na)], condition_na[-1])
                                        start_na <- which(condition_na == "FT") + 1
                                        if (condition_na[1] == "TT" | 
                                            condition_na[1] == "TF"){
                                                    start_na <- c(1, start_na)}
                                        end_na <- which(condition_na == "TF")
                                        if (condition_na[length(condition_na)] == "TT" | 
                                            condition_na[length(condition_na)] == "FT"){
                                                        end_na <- c(end_na, length(condition_na) + 1)}
                                        if (length(start_na) > 0){res <- end_na - start_na + 1} else {res <- NA}
                                } else {res <- NA}
                                return(res)},
                        mc.cores = cores)}
        consecutive_na <- unname(unlist(consecutive_na))
        consecutive_na <- consecutive_na[!is.na(consecutive_na)]
        result <- data.frame(
                "window_size" = max(consecutive_na):1,
                "min_log2ratio_na" = rep(0, times = max(consecutive_na)),
                "max_log2ratio_na" = rep(0, times = max(consecutive_na)))
        consecutive_na <- table(consecutive_na)
        for (w_ind in 1:nrow(result)){
                w <- as.numeric(result[w_ind,1])
                num_w <- unname(consecutive_na[match(w, as.numeric(as.character(names(consecutive_na))))])
                if (!is.na(num_w)){
                        coef_min <- 1:w
                        coef_max <- seq(from = 2, to = w+2, by = 2)
                        coef_max <- c(coef_max, rep(w+2, times = w - length(coef_max)))
                        result[w_ind:nrow(result), 2] <- result[w_ind:nrow(result), 2]  + num_w*coef_min
                        result[w_ind:nrow(result), 3] <- result[w_ind:nrow(result), 2]  + num_w*coef_max}}
        result <- result[nrow(result):1,]
        if (is.na(limit_size)){
                total_log2ratio <- length(data_list) - (2*result[,1] + 1)
                total_log2ratio <- total_log2ratio*(total_log2ratio - 1)/2 + total_log2ratio
        } else {
                if (limit_size + 1 >= length(data_list)){
                        warning(paste0("The parameter limit_size exceeds the number of bins available. Non-restricted matrix is used."))
                        total_log2ratio <- length(data_list) - (2*result[,1] + 1)
                        total_log2ratio <- total_log2ratio*(total_log2ratio - 1)/2 + total_log2ratio
                } else {
                        total_log2ratio <- (length(data_list) - limit_size)*(limit_size - 2*result[,1])}}
        result$min_log2ratio_non_na <- pmax(total_log2ratio - result[,3],0)
        result$max_log2ratio_non_na <- pmax(total_log2ratio - result[,2],0)
        return(result)
}
