
prefix <- "dpnII_raw"
select_limit_size <- 1000

chr_vector <- c("2L", "2R", "3L", "3R", "4", "X", "Y")

############################################ 
###########  PERFORM WINDOW SIZE DIAGNOSTICS
############################################ 

directory <- paste0("~/cortado/", prefix, "/tads/processing/")

############################################ COrTADo start

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, "_", region, "_list_per_row_rds")
        hic_region_list_per_row <- readRDS(paste0(directory, rds_name))
        rm(rds_name)

        rds_name <- paste0("BG3_", prefix, "_", region, "_window_size_diagnostics_table_per_row_rds")
        window_size_diagnostics_table_per_row <- window_size_diagnostics(
                data_list = hic_region_list_per_row, limit_size = select_limit_size, cores = select_n_cores)        
        saveRDS(window_size_diagnostics_table_per_row, paste0(direcotry, rds_name))
        rm(window_size_diagnostics_table_per_row, hic_region_list_per_row, rds_name)
}

############################################ COrTADo end

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, "_", region, "_list_per_col_rds")
        hic_region_list_per_col <- readRDS(paste0(directory, rds_name))
        rm(rds_name)

        rds_name <- paste0("BG3_", prefix, "_", region, "_window_size_diagnostics_table_per_col_rds")
        window_size_diagnostics_table_per_col <- window_size_diagnostics(
                data_list = hic_region_list_per_col, limit_size = select_limit_size, cores = select_n_cores)        
        saveRDS(window_size_diagnostics_table_per_col, paste0(direcotry, rds_name))
        rm(window_size_diagnostics_table_per_col, hic_region_list_per_col, rds_name)
}
