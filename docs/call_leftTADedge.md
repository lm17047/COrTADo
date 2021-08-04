## call_leftTADedge()

### Description



### Usage

```{r}
call_leftTADedge(
    data_list, window_size = 10, prob_limit = 0.001, 
    replace_zero = TRUE, bandwidth_size = 10, quant_limit = 0.90, 
    show_all = TRUE, do_onesided = TRUE, do_correction = TRUE, cores = 10)
```

### Required arguments

`data_list`  
The output of `compute_log2mean()` function produced at the previous stage. A list object where each element is a data.frame indicating log2 mean ratios computed at specified column position and extracted row-wise. 

`window_size`  
A non-negative integer. The desired length of validation test windows.

### Optional arguments

`replace_zero`  
Logical value indicating whether the log2 mean ratio resulting in NA should be replaced with zero.  
Default: TRUE.

`quant_limit`  
A quantile restriction parameter. Removes log2 mean ratios at each column position that are lower than specified quantile value. Ignored if NA.   
Default: NA.

`bandwidth_size`  
A kernel bandwidth smoothing parameter (see `ksmooth()` R documentation for more details). If NA the simple average of log2 mean ratio at the column position is used.   
Default: NA.

`prob_limit`  
P-value threshold for the probability of a column candidate position to be a local maximum.   
Default: 0.001.

`show_all`  
Logical value indicating whether all local maxima candidate positions are reported, not only significant ones.  
Default: TRUE.

`do_onesided`  
Logical value indicating the alternative hypothesis at Wilcoxon Rank-Sum test to validate local maxima. If TRUE, the alternative is "greater", the validation test checks whether the distribution within middle window is greater than within left and right windows. If FALSE, the alternative is "two.sided", the validation test checks whether the distribution within middle window is different from the distributions within left and right windows.  
Default: TRUE.

`do_correction`  
Logical value indicating whether the Bonferonni multiple testing correction is applied.  
Default: TRUE.

`cores`  
A non-negative integer. The number of cores to use to run processes in parallel. If NA, apply non-parallel computing.   
Default: NA.

### Output

A data.frame object containing four columns if `show_all = TRUE` and three columns otherwise: `col_locus` represents the left TAD edge locus name, `col_index` represents the left TAD edge locus index (column index within Hi-C matrix), `p_adj` represents the probability of the column position to be local maximum and `is_significant` represents validated local maximum if `1` and non-validated column candidate position if `0`. Note that the locus names are stored in the following format `chr_start_end`.

