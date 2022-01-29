transform_list2hicmatrix <- function(data_list, bin_first = NA, bin_last = NA, direction){
        if (is.na(bin_first)){bin_first <- 1}
        if (is.na(bin_last)){bin_first <- length(data_list)}
# Create an index for each pair of loci with signal in data_list
        bins_n <- bin_last - bin_first + 1
        locus <- names(data_list)
        if (direction == "row"){
                contact_table <- data.frame(
                        "row_ind" = rep(bin_first:bin_last, times = unname(unlist(lapply(data_list, FUN = nrow)))),
                        "col_ind" = unname(unlist(lapply(data_list, FUN = function(x){x[,2]}))),
                        "contact" = unname(unlist(lapply(data_list, FUN = function(x){x[,3]}))))
                contact_table <- contact_table[!is.na(contact_table[,2]), ]}
        if (direction == "column"){
                contact_table <- data.frame(
                        "row_ind" = unname(unlist(lapply(data_list, FUN = function(x){x[,2]}))),
                        "col_ind" = rep(bin_first:bin_last, times = unname(unlist(lapply(data_list, FUN = nrow)))),
                        "contact" = unname(unlist(lapply(data_list, FUN = function(x){x[,3]}))))
                contact_table <- contact_table[!is.na(contact_table[,1]), ]}
        contact_table$locus <- paste0(contact_table[,1], "_", contact_table[,2])
# Fill matrix with values
        result <- paste0(rep(bin_first:bin_last, each = bins_n), "_", rep(bin_first:bin_last, times = bins_n))        
        result <- contact_table[match(result, contact_table[,4]),3]
        result <- matrix(result, nrow = bins_n, ncol = bins_n, byrow = TRUE)
        colnames(result) <- locus
        rownames(result) <- locus
# Return the resulting matrix
        return(result)
}
