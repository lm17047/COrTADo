## compute_log2mean()

### Description

Compute the measure called log2 mean ratio that is further used to identify the difference between the mean values on the regions that are on the right and on the left from each Hi-C bin. For each column position (where appropriate), we extract the interactions within the window of selected size on the right-hand-side and on the left-hand-side and compute mean values withing these windows. Then, we divide right-hand-side mean over left-hand-side mean and take log2 of the result to get the log2 mean ratio at the analysed column position. 

### Usage

```{r}
compute_log2mean(data_list, replace_zero = FALSE, window_size = NA, cores = NA)
```

### Required arguments

`data_list`  
The output of `transform_hictable2list()` function produced at the previous stage. A list object where each element is a data.frame indicating interactions between pairs of loci extracted row-wise.

`window_size`  
A non-negative integer. The desired length of Moving Average window that is used at log2 mean ratio computation. Has to be specified consistently with `window_size` used in `transform_hictable2list()` function at previous stage. 

### Optional arguments

`replace_zero`  
Logical value indicating whether the log2 mean ratio resulting in NA should be replaced with zero.  
Default: FALSE. 

`cores`  
A non-negative integer. The number of cores to use to run processes in parallel. If NA, apply non-parallel computing.  
Default: NA.

### Output

A list where each element named as row locus and contains a data.frame with three columns: `col_locus` is the column locus name, `col\_ind` is the column locus index within the Hi-C matrix and `log2mean` is the log2 mean ratio computed at column locus position at specified row locus. If the specific data.frame contains only single row (NA, NA, NA) it means that the number of column positions within the row is not enough to compute the log2 mean ratios. Note that the locus names are stored in the following format `chr_start_end`.

