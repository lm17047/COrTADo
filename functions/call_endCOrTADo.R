call_endCOrTADo <- function(data_list, replace_zero = FALSE, window_size, bandwidth_size = NA, do_weighted = TRUE, prob_limit = NA, es_limit = NA, do_onesided = TRUE, test_depth_step, do_prob_correction = FALSE, correction_method = NA, cores = NA){
# Transform log2 mean ratio stored column-wise into row-wise table        
        row_table <- data.frame(
                "row_locus" = unname(unlist(lapply(data_list, "[", 1))),
                "row_ind" = unname(unlist(lapply(data_list, "[", 2))),
                "col_locus" = rep(names(data_list), unname(unlist(lapply(data_list, FUN = nrow)))),
                "col_ind" = rep(1:length(data_list), unname(unlist(lapply(data_list, FUN = nrow)))),
                "log2mean" = unname(unlist(lapply(data_list, "[", 3))))
        row_table[, 1] <-  as.character(row_table[, 1])
        row_table[, 3] <-  as.character(row_table[, 3])
        row_table <- row_table[!is.na(row_table[,1]),]
        row_table <- row_table[order(row_table[,2]),]
# Replace not available log2 mean ratios with zeros or remove 
        if (replace_zero){row_table[is.na(row_table[,5]), 5] <- 0}
# Call local maxima: smooth, define candidate positions and validate through the test
        if (do_weighted){weights <- 1/(row_table[,2] - row_table[,4])} else {weights <- rep(1, times = nrow(row_table))}
        data_smooth <- data.frame(
                        "X" = unique(row_table[,2]),
                        "Y" = unname(unlist(lapply(X = split(weights*row_table[,5], row_table[,2]), FUN = mean, na.rm = TRUE))))
        if (!is.na(bandwidth_size)){
                data_smooth <- .call_smooth(
                        data_table = data.frame("X" = data_smooth[,1], "Y" = data_smooth[,2]),
                        bandwidth_size = bandwidth_size)}
        data_min <- .call_minima(data_table = data_smooth, window_size = window_size)
        if (is.na(cores)){
                data_prob <- lapply(
                        X = data_min[,1], FUN = .valid_test_min, data_table = row_table, 
                        window_size = window_size, do_onesided = do_onesided, test_depth_step = test_depth_step,
                        start_limit = min(row_table[,2]), end_limit = max(row_table[,2]))
        } else {
                data_prob <- parallel::mclapply(
                        data_min[,1], .valid_test_min, data_table = row_table, 
                        window_size = window_size, do_onesided = do_onesided, test_depth_step = test_depth_step,
                        start_limit = min(row_table[,2]), end_limit = max(row_table[,2]))}
# Combine the results                       
        result <- data.frame(
                "row_locus" = names(data_list)[data_min[,1]],
                "row_index" = data_min[,1], 
                "depth" = unname(unlist(data_prob))[c(F,F,T,F)],
                "strength" = unname(unlist(data_prob))[c(F,F,F,T)],
                "effect_size" = unname(unlist(data_prob))[c(F,T,F,F)],
                "p_value" = unname(unlist(data_prob))[c(T,F,F,F)])
        result <- result[!duplicated(result[,2]),]
# Do multi-test correction and return the results                          
        result$p_adj <- NA
        if (do_prob_correction){
                result[,7] <- p.adjust(result[,6], method = "bonferroni")
                p_index <- 7} else {p_index <- 6}
        if (is.na(prob_limit) & is.na(es_limit)){
                result$is_significant <- NA
        } else {
                if (!is.na(prob_limit)){
                        prob_condition <- result[,p_index] <= prob_limit & !is.na(result[,p_index])
                } else {prob_condition <- rep(TRUE, times = nrow(prob_condition))}
                if (!is.na(es_limit)){
                        es_condition <- result[,5] >= es_limit & !is.na(result[,5])
                } else {es_condition <- rep(TRUE, times = nrow(prob_condition))}
                result$is_significant <- 0
                result[prob_condition & es_condition] <- 1}
        return(result)
}

# Call local minima from set of consecutive values
.call_minima <- function(data_table, window_size){
        data_table_na_replaced <- data_table
        data_table_na_replaced[is.na(data_table_na_replaced)] <- 0
        Y1 <- data_table_na_replaced[,2]
        Y1 <- Y1[-c(length(Y1) - 1, length(Y1))]
        Y2 <- data_table_na_replaced[,2]
        Y2 <- Y2[-c(1, length(Y2))]
        Y3 <- data_table_na_replaced[,2]
        Y3 <- Y3[-c(1, 2)]
        position <- which(Y1 > Y2 & Y2 < Y3) + 1
        position <- position[!is.na(data_table[position,2])]
        select_condition <- unname(unlist(lapply(X = position, FUN = function(x){
                window_start <- max(1, x - window_size)
                window_end <- min(x + window_size, nrow(data_table))
                window_min <- c(window_start:window_end)[
                        which.min(data_table[window_start:window_end,2])]
                res <- (x == window_min)
                return(res)})))
        position <- position[select_condition]
        result <- data_table[position,]
        return(result)
}

# Validation test at x position
.valid_test_min <- function(x, window_size, data_table, start_limit, end_limit, do_onesided, test_depth_step){
        if (do_onesided){
                alternative <- "less"
                factor <- 1
        } else {
                alternative <- "two.sided"
                factor <- 2}
        start_middle <- x - (window_size %/% 2)
        end_middle <- x + (window_size %/% 2)
        start_left <- start_middle - ((window_size %/% 2)*2 + 1)
        end_left <- start_middle - 1
        start_right <- end_middle + 1
        end_right <- end_middle + ((window_size %/% 2)*2 + 1)
        if (start_left >= start_limit & end_right <= end_limit){
                top_limit <- max(data_table[data_table[,2] == start_left, 4])
                bottom_limit  <- min(data_table[data_table[,2] == end_right, 4])
                if (!is.na(top_limit) & !is.na(bottom_limit)){
                        data_matrix <- matrix(
                                data_table[
                                        data_table[,2] >= start_left & data_table[,2] <= end_right &
                                        data_table[,4] <= top_limit & data_table[,4] >= bottom_limit, 5],
                                nrow = top_limit - bottom_limit + 1, ncol = end_right - start_left + 1, byrow = FALSE)
                        data_matrix <- data_matrix[nrow(data_matrix):1,]
                        if (nrow(data_matrix) >= test_depth_step){
                                class <- rep(1:(nrow(data_matrix)%/% test_depth_step), each = test_depth_step)
                                data_left_list <- split(data_matrix[1:length(class), 1:(ncol(data_matrix)/3)], class)
                                data_middle_list <- split(data_matrix[1:length(class), (ncol(data_matrix)/3+1):(2*ncol(data_matrix)/3)], class)
                                data_right_list <- split(data_matrix[1:length(class), (2*ncol(data_matrix)/3+1):ncol(data_matrix)], class)
                        } else {
                                data_left_list <- list(data_matrix[, 1:(ncol(data_matrix)/3)])
                                data_middle_list <- list(data_matrix[, (ncol(data_matrix)/3+1):(2*ncol(data_matrix)/3)])
                                data_right_list <- list(data_matrix[, (2*ncol(data_matrix)/3+1):ncol(data_matrix)])}
                        result <- unname(unlist(lapply(X = 1:length(data_middle_list), FUN = function(test_ind){
                                data_set_left <- unname(unlist(data_left_list[1:test_ind]))
                                data_set_middle <- unname(unlist(data_middle_list[1:test_ind]))
                                data_set_right <- unname(unlist(data_right_list[1:test_ind]))
                                n_left <- length(data_set_left[!is.na(data_set_left)])
                                n_middle <- length(data_set_middle[!is.na(data_set_middle)])
                                n_rigth <- length(data_set_right[!is.na(data_set_right)])
                                if (n_left > 0 & n_middle > 0 & n_rigth > 0){
                                        p_left <- wilcox.test(
                                                x = data_set_middle, y = data_set_left, 
                                                alternative = alternative, exact = FALSE)$p.value
                                        es_left <- abs(qnorm(p_left/factor))/sqrt(n_middle + n_left)
                                        p_right <- wilcox.test(
                                                x = data_set_middle, y = data_set_right,
                                                alternative = alternative, exact = FALSE)$p.value
                                        es_right <- abs(qnorm(p_right/factor))/sqrt(n_middle + n_rigth)
                                        res <- c( 
                                                min(p_left, p_right),
                                                c(es_left, es_right)[which.min(c(p_left, p_right))])                                        
                                } else {res <- c(NA,NA)}
                                return(res)})))
                        if (all(is.na(result))){
                                result <- c(NA, NA, NA, NA)
                        } else {
                                p <- unname(unlist(result))[c(T,F)]
                                es <- unname(unlist(result))[c(F,T)]
                                if (length(p) > 1){
                                        es_fall <- es[-length(es)] - es[-1]
                                        es_fall[1:which.max(es)] <- NA
                                        if (all(is.na(es_fall))){select_index <- which.max(es)} else {select_index <- which.max(es_fall)}
                                } else {select_index <- 1}
                                if (p[select_index] == Inf | p[select_index] == -Inf | es[select_index] == Inf | es[select_index] == -Inf){
                                        result <- c(NA, NA, NA, NA)
                                } else {
                                        result <- c(
                                                p[select_index], es[select_index],
                                                x - top_limit + select_index*test_depth_step,
                                                mean(unname(unlist(data_middle_list[1:select_index])), na.rm = TRUE))}}
                } else {result <- c(NA, NA, NA, NA)}
        } else {result <- c(NA, NA, NA, NA)}
        return(result)
}
