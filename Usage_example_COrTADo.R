
############################################ 
############################  SET PARAMETERS
############################################ 

prefix <- "dpnII_raw"
input_tsv <- paste0("~/cortado/", prefix, "/BG3_WT_merged_hic_matrix_dpnII.tsv")

select_n_cores <- 30
select_limit_size <- 1000
select_window_size <- 10

select_correction <- "fdr"
select_pval <- 0.01
select_es_weak1 <- 0.3
select_stength_weak1 <- 0.2
select_es_weak2 <- 0.1
select_stength_weak2 <- 0.4
select_es_strong1 <- 0.5
select_stength_strong1 <- 0.4
select_es_strong2 <- 0.3
select_stength_strong2 <- 0.6

############################################ 
######################  REGIONS & TIME TABLE
############################################ 

chr_vector <- c("2L", "2R", "3L", "3R", "4", "X", "Y")

timing_names <- c(
        "data_extraction", 
        "data_transformation_per_row", "data_transformation_per_col", 
        "compute_log2mean_per_row", "compute_log2mean_per_col",
        "call_startCOrTADo", "call_endCOrTADo")
timing_table <- matrix(
        rep(NA, times = length(timing_names)*length(chr_vector)),
        ncol = length(timing_names), nrow = length(chr_vector))
timing_table <- as.data.frame(timing_table)
colnames(timing_table) <- timing_names
rownames(timing_table) <- chr_vector

############################################ 
#######################  TRANSFORM INTO LIST
############################################ 

directory <- paste0("~/cortado/", prefix, "/tads/processing/")

df <- read.table(file = input_tsv, header = FALSE, sep = "\t")

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        start_time <- Sys.time()
        if (!is.na(region)){
                chr <- as.character(unlist(strsplit(region, ":"))[1])
                start <- as.numeric(unlist(strsplit(unlist(strsplit(region, ":"))[2], "-"))[1])
                end <- as.numeric(unlist(strsplit(unlist(strsplit(region, ":"))[2], "-"))[2])
                if (is.na(start) | is.na(end)){
                        warning(paste0("Start or end are not provided, whole chromosome chr", region, " is processed."))
                        select_df <- df[df[,1] == chr & df[,4] == chr, ]
                } else {
                        select_df <- df[
                                df[,1] == chr & df[,4] == chr &
                                df[,2] >= start & df[,5] >= start &
                                df[,3] <= end & df[,6] <= end, ]}}
        end_time <- Sys.time()
        timing_table[chr_index, 1] <- as.numeric(difftime(end_time, start_time, units= "mins"))
        rm(start_time, end_time)
        
############################################ COrTADo start

        rds_name <- paste0("BG3_", prefix, region, "_list_per_row_rds")
        start_time <- Sys.time()
        hic_region_list_per_row <- transform_hictable2list(
                data_table = select_df, direction = "row", replace_zero = FALSE, resolution = NA,
                limit_size = select_limit_size, cores = select_n_cores)
        end_time <- Sys.time()
        timing_table[chr_index, 2] <- as.numeric(difftime(end_time, start_time, units= "mins"))
        saveRDS(hic_region_list_per_row, rds_name)
        rm(hic_region_list_per_row, rds_name, start_time, end_time)
        
############################################ COrTADo end

        rds_name <-  paste0("BG3_", region, "_list_per_col_rds")
        start_time <- Sys.time()
        hic_region_list_per_col <- transform_hictable2list(
                data_table = select_df, direction = "column", replace_zero = FALSE, resolution = NA,
                limit_size = select_limit_size, cores = select_n_cores)
        end_time <- Sys.time()
        timing_table[chr_index, 3] <- as.numeric(difftime(end_time, start_time, units= "mins"))
        saveRDS(hic_region_list_per_col, rds_name)
        rm(hic_region_list_per_col, rds_name, start_time, end_time)

        rm(select_df, start_time, end_time)
}

rm(df, start_time, end_time)
saveRDS(timing_table, "timing_table_rds")

############################################ END OF SECTION



############################################ 
#################  COMPUTE LOG2MEAN RATIOS
############################################ 

directory <- paste0("~/cortado/", prefix, "/tads/processing/")

############################################ COrTADo start

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, region, "_list_per_row_rds")
        hic_region_list_per_row <- readRDS(paste0(directory, rds_name))
        rm(rds_name)

        rds_name <- paste0("BG3_", prefix, region, "_log2mean_list_per_row_rds")
        start_time <- Sys.time()
        log2mean_list_per_row <- compute_log2mean(
                data_list = hic_region_list_per_row, window_size = select_window_size, 
                replace_zero = FALSE, cores = select_n_cores)        
        end_time <- Sys.time()
        timing_table[chr_index,4] <- as.numeric(difftime(end_time, start_time, units= "mins"))
        saveRDS(log2mean_list_per_row, paste0(rds_name))
        rm(log2mean_list_per_row, hic_region_list_per_row, rds_name, start_time, end_time)
}

############################################ COrTADo end

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, region, "_list_per_col_rds")
        hic_region_list_per_col <- readRDS(paste0(directory, rds_name))
        rm(rds_name)

        rds_name <- paste0("BG3_", prefix, region, "_log2mean_list_per_col_rds")
        start_time <- Sys.time()
        log2mean_list_per_col <- compute_log2mean(
                data_list = hic_region_list_per_col, window_size = select_window_size, 
                replace_zero = FALSE, cores = select_n_cores)        
        end_time <- Sys.time()
        timing_table[chr_index, 5] <- as.numeric(difftime(end_time, start_time, units= "mins"))
        saveRDS(log2mean_list_per_row, paste0(rds_name))
        rm(log2mean_list_per_col, hic_region_list_per_col, rds_name, start_time, end_time)
}

saveRDS(timing_table, "timing_table_rds")

############################################ END OF SECTION



############################################ 
########################  CALL START AND END 
############################################ 

directory <- paste0("~/cortado/", prefix, "/tads/processing/")

############################################ COrTADo start

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, region, "_log2mean_list_per_row_rds")
        log2mean_list_per_row <- readRDS(paste0(directory, rds_name))
        rm(rds_name)

        rds_name <- paste0("BG3_", prefix, region, "_startCOrTADo_table_rds")
        start_time <- Sys.time()
        startCOrTADo_table <- call_startCOrTADo(
                data_list = log2mean_list_per_row, replace_zero = FALSE, 
                window_size = select_window_size, bandwidth_size = NA, 
                do_weighted = TRUE, prob_limit = NA, es_limit = NA, 
                do_onesided = TRUE, do_prob_correction = FALSE, correction_method = NA,
                test_depth_step = select_window_size, cores = select_n_cores)
        end_time <- Sys.time()
        timing_table[chr_index, 6] <- as.numeric(difftime(end_time, start_time, units= "mins"))
        saveRDS(startCOrTADo_table, paste0(rds_name))
        rm(log2mean_list_per_row, startCOrTADo_table, rds_name, start_time, end_time)
}

############################################ COrTADo end

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, region, "_log2mean_list_per_col_rds")
        log2mean_list_per_row <- readRDS(paste0(directory, rds_name))
        rm(rds_name)

        rds_name <- paste0("BG3_", prefix, region, "_endCOrTADo_table_rds")
        start_time <- Sys.time()
        endCOrTADo_table <- call_endCOrTADo(
                data_list = log2mean_list_per_col, replace_zero = FALSE, 
                window_size = select_window_size, bandwidth_size = NA, 
                do_weighted = TRUE, prob_limit = NA, es_limit = NA, 
                do_onesided = TRUE, do_prob_correction = FALSE, correction_method = NA,
                test_depth_step = select_window_size, cores = select_n_cores)
        end_time <- Sys.time()
        timing_table[chr_index, 7] <- as.numeric(difftime(end_time, start_time, units= "mins"))
        saveRDS(endCOrTADo_table, paste0(rds_name))
        rm(log2mean_list_per_col, endCOrTADo_table, rds_name, start_time, end_time)
}

saveRDS(timing_table, "timing_table_rds")

############################################ END OF SECTION



############################################ 
#########  AGGREGATE ACROSS CHROMOSOMES
############################################ 

directory <- paste0("~/cortado/", prefix, "/tads/processing/")

############################################ COrTADo start

chr_length <- 0
startTADedge_table <- c(NA, NA, NA, NA, NA)

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, "_", region, "_startCOrTADo_table_rds")
        startCOrTADo_region <- readRDS(paste0(directory, rds_name))
        startCOrTADo_region[,2] <- startCOrTADo_region[,2] + chr_length
        startTADedge_table <- rbind(startTADedge_table, startCOrTADo_region)
        rm(startCOrTADo_region, rds_name)
        
        rds_name <-  paste0("BG3_", prefix, "_", region,  "_log2mean_list_per_row_rds")
        chr_length <- chr_length + length(readRDS(paste0(directory, rds_name)))
        rm(rds_name)
}

startTADedge_table <- startTADedge_table[-1,]
rds_name <- paste0(directory, "BG3_", prefix, "_startCOrTADo_table_combined_rds")
saveRDS(startTADedge_table, rds_name)
rm(startTADedge_table, chr_length, rds_name)

############################################ COrTADo end

chr_length <- 0
endTADedge_table <- c(NA, NA, NA, NA, NA)

for (chr_index in 1:length(chr_vector)){
region <- chr_vector[chr_index]
print(paste0("Region ", region, " is processing."))
        rds_name <- paste0("BG3_", prefix, "_", region, "_endCOrTADo_table_rds")
        endCOrTADo_region <- readRDS(paste0(directory, rds_name))
        endCOrTADo_region[,2] <- endCOrTADo_region[,2] + chr_length
        endTADedge_table <- rbind(endTADedge_table, endCOrTADo_region)
        rm(rds_name, endCOrTADo_region)
        
        rds_name <-  paste0("BG3_", prefix, "_", region,  "_log2mean_list_per_col_rds")
        chr_length <- chr_length + length(readRDS(paste0(directory, rds_name)))
        rm(rds_name)
}

endTADedge_table <- endTADedge_table[-1,]
rds_name <- paste0(directory, "BG3_", prefix, "_endCOrTADo_table_combined_rds")
saveRDS(endTADedge_table, rds_name)
rm(endTADedge_table, chr_length, rds_name)

############################################ END OF SECTION



############################################ 
#############  ALL, WEAK & STRONG SELECTION
############################################ 

directory <- paste0("~/cortado/", prefix, "/tads/processing/")

############################################ COrTADo start

rds_name <- paste0("BG3_", prefix, "_startCOrTADo_table_combined_rds")
startCOrTADo_combined <- readRDS(paste0(directory, rds_name))
startCOrTADo_combined[,7] <- p.adjust(startCOrTADo_combined[,6], method = select_correction)
startCOrTADo_combined <- startCOrTADo_combined[
        !is.na(startCOrTADo_combined[,3]) &
        !is.na(startCOrTADo_combined[,4]) & startCOrTADo_combined[,4]>= 0 &
        !is.na(startCOrTADo_combined[,5]) & 
        !is.na(startCOrTADo_combined[,7]),]
rm(rds_name)

startCOrTADo_combined[,8] <- 0
startCOrTADo_combined[startCOrTADo_combined[,7] <= select_pval,8] <- 1
startCOrTADo_combined_all <- startCOrTADo_combined[startCOrTADo_combined[,8] == 1,-8]
rds_name <- paste0(directory, "BG3_", prefix, "_startCOrTADo_table_combined_all_rds")
saveRDS(startCOrTADo_combined_all, rds_name)
rm(startCOrTADo_combined_all, rds_name)

startCOrTADo_combined[,8] <- 0
startCOrTADo_combined[
        startCOrTADo_combined[,7] <= select_pval &
        ((startCOrTADo_combined[,4] >= select_stength_weak1 & startCOrTADo_combined[,5] >= select_es_weak1)|
        (startCOrTADo_combined[,4] >= select_stength_weak2 & startCOrTADo_combined[,5] >= select_es_weak2)),8] <- 1
startCOrTADo_combined_weak <- startCOrTADo_combined[startCOrTADo_combined[,8] == 1,-8]
rds_name <- paste0(directory, "BG3_", prefix, "_startCOrTADo_table_combined_weak_rds")
saveRDS(startCOrTADo_combined_weak, rds_name)
rm(startCOrTADo_combined_weak, rds_name)

startCOrTADo_combined[,8] <- 0
startCOrTADo_combined[
        startCOrTADo_combined[,7] <= select_pval &
        ((startCOrTADo_combined[,4] >= select_stength_strong1 & startCOrTADo_combined[,5] >= select_es_strong1)|
        (startCOrTADo_combined[,4] >= select_stength_strong2 & startCOrTADo_combined[,5] >= select_es_strong2)),8] <- 1
startCOrTADo_combined_strong <- startCOrTADo_combined[startCOrTADo_combined[,8] == 1,-8]
rds_name <- paste0(directory, "BG3_", prefix, "_startCOrTADo_table_combined_strong_rds")
saveRDS(startCOrTADo_combined_strong, rds_name)
rm(startCOrTADo_combined_strong, rds_name)

############################################ COrTADo end

rds_name <- paste0("BG3_", prefix, "_endCOrTADo_table_combined_rds")
endCOrTADo_combined <- readRDS(paste0(directory, rds_name))
endCOrTADo_combined[,7] <- p.adjust(endCOrTADo_combined[,6], method = select_correction)
endCOrTADo_combined <- endCOrTADo_combined[
        !is.na(endCOrTADo_combined[,3]) &
        !is.na(endCOrTADo_combined[,4]) & endCOrTADo_combined[,4]<= 0 &
        !is.na(endCOrTADo_combined[,5]) & 
        !is.na(endCOrTADo_combined[,7]),]
endCOrTADo_combined[,4] <- abs(endCOrTADo_combined[,4])
rm(rds_name)

endCOrTADo_combined[,8] <- 0
endCOrTADo_combined[endCOrTADo_combined[,7] <= select_pval,8] <- 1
endCOrTADo_combined_all <- endCOrTADo_combined[endCOrTADo_combined[,8] == 1,-8]
rds_name <- paste0(directory, "BG3_", prefix, "_endCOrTADo_table_combined_all_rds")
saveRDS(endCOrTADo_combined_all, rds_name)
rm(endCOrTADo_combined_all, rds_name)

endCOrTADo_combined[,8] <- 0
endCOrTADo_combined[
        endCOrTADo_combined[,7] <= select_pval &
        ((endCOrTADo_combined[,4] >= select_stength_weak1 & endCOrTADo_combined[,5] >= select_es_weak1)|
        (endCOrTADo_combined[,4] >= select_stength_weak2 & endCOrTADo_combined[,5] >= select_es_weak2)),8] <- 1
endCOrTADo_combined_weak <- endCOrTADo_combined[endCOrTADo_combined[,8] == 1,-8]
rds_name <- paste0(directory, "BG3_", prefix, "_endCOrTADo_table_combined_weak_rds")
saveRDS(endCOrTADo_combined_weak, rds_name)
rm(endCOrTADo_combined_weak, rds_name)

endCOrTADo_combined[,8] <- 0
endCOrTADo_combined[
        endCOrTADo_combined[,7] <= select_pval &
        ((endCOrTADo_combined[,4] >= select_stength_strong1 & endCOrTADo_combined[,5] >= select_es_strong1)|
        (endCOrTADo_combined[,4] >= select_stength_strong2 & endCOrTADo_combined[,5] >= select_es_strong2)),8] <- 1
endCOrTADo_combined_strong <- endCOrTADo_combined[endCOrTADo_combined[,8] == 1,-8]
rds_name <- paste0(directory, "BG3_", prefix, "_endCOrTADo_table_combined_strong_rds")
saveRDS(endCOrTADo_combined_strong, rds_name)
rm(endCOrTADo_combined_strong, rds_name)

rm(directory, startCOrTADo_combined, endCOrTADo_combined)

############################################ END OF SECTION
