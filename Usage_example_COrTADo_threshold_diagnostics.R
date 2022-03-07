
prefix <- "dpnII_raw"

setwd("~/cortado/")
library(png)

############################################ Load data

input_directory <- paste0("~/cortado/", prefix, "/tads/processing/")

input_file <- paste0("BG3_", prefix, "_startCOrTADo_table_combined_rds")
data_start_full <- readRDS(paste0(input_directory, input_file))
data_start_full$p_adj <- p.adjust(data_start_full$p_value, method = select_correction)
data_start_full <- data_start_full[
        !is.na(data_start_full$strength) & data_start_full$strength >= 0 &
        !is.na(data_start_full$effect_size) & 
        !is.na(data_start_full$p_adj), ]
data_start_full <- data_start_full[data_start_full$p_adj <= 0.01,]

input_file <- paste0("BG3_", prefix, "_endCOrTADo_table_combined_rds")
data_end_full <- readRDS(paste0(input_directory, input_file))
data_end_full$p_adj <- p.adjust(data_end_full$p_value, method = select_correction)
data_end_full <- data_end_full[
        !is.na(data_end_full$strength) & data_end_full$strength <= 0 &
        !is.na(data_end_full$effect_size) & 
        !is.na(data_end_full$p_adj), ]
data_end_full <- data_end_full[data_end_full$p_adj <= 0.01,]
data_end_full[,4] <- abs(data_end_full$strength)

############################################ Strength vs effect size

pdf(paste0("~/cortado/", prefix, "/tads/analysis/BG3_", prefix, "_startCOrTADo_threshold_diagnostics_volcano_es_vs_strength.pdf"), width = 5, height = 6)
        par(mar = c(6.0, 6.0, 4.0, 4.0))
        par(mgp = c(0.5, 2.0, 1.0))
        x_breaks <- c(0,0.2,0.4,0.6,0.8,1.0,1.2)
        y_breaks <- c(0,0.1,0.3,0.5,1.0)
        plot(
                x = 0, y = 0, type = "n", xlab = "", ylab = "", axes = FALSE, yaxs="i", xaxs="i",
                xlim = c(x_breaks[1],x_breaks[length(x_breaks)]), 
                ylim = c(y_breaks[1],y_breaks[length(y_breaks)]))
        axis(side = 1, at = x_breaks, labels = x_breaks)
        axis(side = 2, at = y_breaks, labels = y_breaks, las = 2)
        coords <- par("usr")
        gx <- grconvertX(coords[1:2], "user", "in")
        gy <- grconvertY(coords[3:4], "user", "in")
        width <- max(gx) - min(gx)
        height <- max(gy) - min(gy)
        tmp <- tempfile()
                png(tmp, width = width, height = height, units = "in", res = 500, bg = "transparent")
                        par(mar = c(0, 0, 0, 0))
                        plot(
                                x = 0, y = 0, type = "n", xlab = "", ylab = "", axes = FALSE, yaxs="i", xaxs="i",
                                xlim = c(x_breaks[1],x_breaks[length(x_breaks)]), 
                                ylim = c(y_breaks[1],y_breaks[length(y_breaks)]))
                        points(
                                x = data_start_full$strength, y = data_start_full$effect_size, 
                                col = "grey", cex = 0.5, pch = 16)
        dev.off()
        panel <- readPNG(tmp)
        rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
                
        x_breaks[length(x_breaks)] <- Inf
        for (i in 1:(length(x_breaks)-1)){
                for (j in 1:(length(y_breaks)-1)){
                       count <- sum(
                            data_start_full$strength >= x_breaks[i] & 
                            data_start_full$strength < x_breaks[i+1] &
                            data_start_full$effect_size >= y_breaks[j] & 
                            data_start_full$effect_size < y_breaks[j+1])
                        perc <- paste0(round(100*count/nrow(data_start_full), 2), "%")
                        if (x_breaks[i+1] != Inf){
                                text(
                                        x = (x_breaks[i]+x_breaks[i+1])/2, 
                                        y = (y_breaks[j]+y_breaks[j+1])/2, 
                                        labels = paste0(count,"\n",perc), cex = 0.5)
                        } else {                                
                                text(
                                        x = (x_breaks[i]+1.2)/2, 
                                        y = (y_breaks[j]+y_breaks[j+1])/2, 
                                        labels = paste0(count,"\n",perc), cex = 0.5)}}}

        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NA)   
        abline(h = c(0.5,0.3,0.1), lwd = 1, lty=2, col = "black")  
        abline(v = c(0.2,0.4,0.6,0.8,1.0,1.2), lwd = 1, lty=2, col = "black")  
dev.off()


pdf(paste0("~/cortado/", prefix, "/tads/analysis/BG3_", prefix, "_endCOrTADo_threshold_diagnostics_volcano_es_vs_strength.pdf"), width = 5, height = 6)
        par(mar = c(6.0, 6.0, 4.0, 4.0))
        par(mgp = c(0.5, 2.0, 1.0))
        x_breaks <- c(0,0.2,0.4,0.6,0.8,1.0,1.2)
        y_breaks <- c(0,0.1,0.3,0.5,1.0)
        plot(
                x = 0, y = 0, type = "n", xlab = "", ylab = "", axes = FALSE, yaxs="i", xaxs="i",
                xlim = c(x_breaks[1],x_breaks[length(x_breaks)]), 
                ylim = c(y_breaks[1],y_breaks[length(y_breaks)]))
        axis(side = 1, at = x_breaks, labels = x_breaks)
        axis(side = 2, at = y_breaks, labels = y_breaks, las = 2)
        coords <- par("usr")
        gx <- grconvertX(coords[1:2], "user", "in")
        gy <- grconvertY(coords[3:4], "user", "in")
        width <- max(gx) - min(gx)
        height <- max(gy) - min(gy)
        tmp <- tempfile()
                png(tmp, width = width, height = height, units = "in", res = 500, bg = "transparent")
                        par(mar = c(0, 0, 0, 0))
                        plot(
                                x = 0, y = 0, type = "n", xlab = "", ylab = "", axes = FALSE, yaxs="i", xaxs="i",
                                xlim = c(x_breaks[1],x_breaks[length(x_breaks)]), 
                                ylim = c(y_breaks[1],y_breaks[length(y_breaks)]))
                        points(
                                x = data_end_full$strength, y = data_end_full$effect_size, 
                                col = "grey", cex = 0.5, pch = 16)
        dev.off()
        panel <- readPNG(tmp)
        rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
                
        x_breaks[length(x_breaks)] <- Inf
        for (i in 1:(length(x_breaks)-1)){
                for (j in 1:(length(y_breaks)-1)){
                       count <- sum(
                            data_end_full$strength >= x_breaks[i] & 
                            data_end_full$strength < x_breaks[i+1] &
                            data_end_full$effect_size >= y_breaks[j] & 
                            data_end_full$effect_size < y_breaks[j+1])
                        perc <- paste0(round(100*count/nrow(data_end_full), 2), "%")
                        if (x_breaks[i+1] != Inf){
                                text(
                                        x = (x_breaks[i]+x_breaks[i+1])/2, 
                                        y = (y_breaks[j]+y_breaks[j+1])/2, 
                                        labels = paste0(count,"\n",perc), cex = 0.5)
                        } else {                                
                                text(
                                        x = (x_breaks[i]+1.2)/2, 
                                        y = (y_breaks[j]+y_breaks[j+1])/2, 
                                        labels = paste0(count,"\n",perc), cex = 0.5)}}}

        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NA)   
        abline(h = c(0.5,0.3,0.1), lwd = 1, lty=2, col = "black")  
        abline(v = c(0.2,0.4,0.6,0.8,1.0,1.2), lwd = 1, lty=2, col = "black")  
dev.off()


############################################ Weak and Strong

count <- sum(
        (data_start_full$strength >= select_stength_weak1 & data_start_full$effect_size >= select_es_weak1)|
        (data_start_full$strength >= select_stength_weak2 & data_start_full$effect_size >= select_es_weak2))
perc <- count/nrow(data_start_full)
print(count)
print(perc)

count <- sum(
        (data_start_full$strength >= select_stength_strong1 & data_start_full$effect_size >= select_es_strong1)|
        (data_start_full$strength >= select_stength_strong2 & data_start_full$effect_size >= select_es_strong2))
perc <- count/nrow(data_start_full)
print(count)
print(perc)

count <- sum(
        (data_end_full$strength >= select_stength_weak1 & data_end_full$effect_size >= select_es_weak1)|
        (data_end_full$strength >= select_stength_weak2 & data_end_full$effect_size >= select_es_weak2))
perc <- count/nrow(data_end_full)
print(count)
print(perc)

count <- sum(
        (data_end_full$strength >= select_stength_strong1 & data_end_full$effect_size >= select_es_strong1)|
        (data_end_full$strength >= select_stength_strong2 & data_end_full$effect_size >= select_es_strong2))
perc <- count/nrow(data_end_full)
print(count)
print(perc)
