
call_leftTADedge <- function(data_list, replace_zero = TRUE, window_size = 10, bandwidth_size = 10, quant_limit = 0.90, show_all = TRUE, prob_limit = 0.001, do_onesided = TRUE, do_correction = TRUE, cores = NA){

# Transform log2 mean ratio stored row-wise into column-wise table

        col_table <- data.frame(
                "col_locus" = unname(unlist(lapply(data_list, "[", 1))),
                "col_ind" = unname(unlist(lapply(data_list, "[", 2))),
                "row_locus" = rep(names(data_list), unname(unlist(lapply(data_list, FUN = nrow)))),
                "row_ind" = rep(1:length(data_list), unname(unlist(lapply(data_list, FUN = nrow)))),
                "log2mean" = unname(unlist(lapply(data_list, "[", 3))))
        col_table[, 1] <-  as.character(col_table[, 1])
        col_table[, 3] <-  as.character(col_table[, 3])
        col_table <- col_table[!is.na(col_table[,1]),]
        col_table <- col_table[order(col_table[,2]),]
        
# Replace not available log2 mean ratios with zeros or remove 

        if (replace_zero){col_table[is.na(col_table[,5]), 5] <- 0}
        if (!replace_zero){col_table <- col_table[!is.na(col_table[,5]), ] }
        
# Leave only values that are higher than quant_limit 

        if (!is.na(quant_limit)){
                col_table[,6] <- rep(
                        unname(unlist(lapply(
                                X = split(col_table[,5], col_table[,2]), 
                                FUN = quantile, prob = quant_limit, na.rm = TRUE))), 
                        unname(unlist(lapply(
                                X = split(col_table[,4], col_table[,2]), 
                                FUN = length))))
                col_table <- col_table[col_table[,5] >= col_table[,6], -6]}
                
# Transform table into list, each element is column-wise log2 mean ratios

        col_list <- split(col_table[,-c(1,2)], col_table[,1])
        col_list <- col_list[match(names(data_list), names(col_list))]
        names(col_list) <- names(data_list)
        empty_cond <- unname(unlist(lapply(col_list, FUN = is.null)))
        col_list[empty_cond] <- lapply(
                X = col_list[empty_cond], 
                FUN = function(x){x <- data.frame("row_locus" = NA, "row_ind" = NA, "log2mean" = NA)}) 
                
# Call local maxima: smooth, define candidate positions, validate through the test and do multi-test correction

        data_smooth <- .call_smooth(
                data_table = data.frame("X" = col_table[,2], "Y" = col_table[,5]), 
                bandwidth_size = bandwidth_size)
        data_max <- .call_maxima(data_table = data_smooth)
        if (is.na(cores)){
                data_prob <- unlist(lapply(
                        X = data_max[,1], FUN = .valid_test, data_list = col_list,
                        window_size = window_size, do_onesided = do_onesided,
                        start_limit = min(col_table[,2]), end_limit = max(col_table[,2])))
        } else {
                data_prob <- parallel::mclapply(
                        data_max[,1], .valid_test, data_list = col_list, 
                        window_size = window_size, do_onesided = do_onesided,
                        start_limit = min(col_table[,2]), end_limit = max(col_table[,2]))
                data_prob <- unlist(data_prob)}
        if (do_correction){data_prob <- p.adjust(data_prob, method = "bonferroni")}
        
# Combine and return results

        result <- data.frame(
                "col_locus" = names(col_list)[data_max[,1]],
                "col_index" = data_max[,1],
                "strength" = data_max[,2],
                "p_adj" = data_prob)
        if (!is.na(prob_limit)){
                result$is_significant <- 0
                result[result[,4] <= prob_limit & !is.na(result[,4]), 5] <- 1}
        if (!show_all){result <- result[result[,5] == 1,-5]}
        return(result)
}

### Supporting functions .call_smooth(), .call_maxima() and .valid_test()

# Call smoothed values at x positions 

.call_smooth <- function(data_table, bandwidth_size){
        result <- ksmooth(
                x = data_table$X, y = data_table$Y, 
                kernel = "normal", bandwidth = bandwidth_size, 
                n.points  = length(unique(data_table$X)))
        result <- data.frame(X = result[[1]], Y = result[[2]])
        return(result)
}

# Call local maxima from set of consecutive values

.call_maxima <- function(data_table){
        Y1 <- data_table[,2]
        Y1 <- Y1[-c(length(Y1) - 1, length(Y1))]
        Y2 <- data_table[,2]
        Y2 <- Y2[-c(1, length(Y2))]
        Y3 <- data_table[,2]
        Y3 <- Y3[-c(1, 2)]
        position <- which(Y1 < Y2 & Y2 > Y3) + 1
        result <- data_table[position,]
        return(result)
}

# Validation test at x position

.valid_test <- function(x, window_size, data_list, start_limit, end_limit, do_onesided){
        start_middle <- ceiling(x - 0.5*window_size)
        end_middle <- floor(x + 0.5*window_size)
        start_left <- start_middle - window_size
        end_left <- start_middle - 1
        start_right <- end_middle + 1
        end_right <- end_middle + window_size
        if (start_left >= start_limit & end_right <= end_limit){
                data_set_middle <- unname(unlist(lapply(
                        data_list[start_middle:end_middle], "[", 3)))
                data_set_left <- unname(unlist(lapply(
                        data_list[start_left:end_left], "[", 3)))
                data_set_right <- unname(unlist(lapply(
                        data_list[start_right:end_right], "[", 3)))
                data_set_middle <- data_set_middle[!is.na(data_set_middle)]
                data_set_left <- data_set_left[!is.na(data_set_left)]
                data_set_right <- data_set_right[!is.na(data_set_right)]
                if (length(data_set_middle) > 0 & length(data_set_left) > 0 & length(data_set_right)){
                        if (do_onesided){alternative <- "greater"} else {alternative <- "two.sided"}
                        p_left <- wilcox.test(
                                x = data_set_middle, y = data_set_left, 
                                alternative = alternative, exact = FALSE)$p.value
                        p_right <- wilcox.test(
                                x = data_set_middle, y = data_set_right,
                                alternative = alternative, exact = FALSE)$p.value
                        result <- max(c(p_left, p_right))
                } else {result <- NA}
        } else {result <- NA}
        return(result)
}

