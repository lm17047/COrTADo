## compute_log2mean()

### Description

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
