## transform_hictable2list()

### Description

Transforms the input Hi-C table into a list where each element represents the interactions happen within single Hi-C row. The Hi-C table contains seven columns where first three represents the interacting row locus (chromosome, start and end positions), next three represents the interacting column locus (chromosome, start and end positions) and the last one represents the interaction frequency (discrete if Hi-C data is non-normalised and continuous if Hi-C data is normalised). The table omits non-interacting pairs, i.e. if reads of fragments, which do or do not interact, are not mapped to the genome, the loci pair is not represented in the table.

In case when we call TADs on a large scale genomic region (chromosome or genome-wise), we do not need the whole row to be extracted as it requires more memory and more processing time - we need only the part of the matrix that is close to the main diagonal and the size of this part is related to the maximum TAD width that we expect to detect. The large maximum TAD width leads to the large regions at the start and the end of chromatin, where TADs cannot be detected. However, when we work with relatively short part of chromosome, like the testing genomic region, it can become critical, as we extermely reduce the number of observations extracted, so we prefer to work with the whole unrestricted lower triangle matrix instead.

<p align="center">
<img src="https://github.com/lm17047/TADedge_calling/blob/active/docs/img/img2.png" width="800">
</p> 

### Usage

```{r}
transform_hictable2list(data_table, replace_zero = FALSE, limit_size = NA, window_size = NA, cores = NA)
```

### Required arguments

`data_table`  
A data.frame object containing seven columns where first three represents the interacting row locus (chromosome, start and end positions), next three represents the interacting column locus (chromosome, start and end positions) and the last one represents the interaction frequency (discrete if Hi-C data is non-normalised and continuous if Hi-C data is normalised).

### Optional arguments

`replace_zero`  
Logical value indicating whether the interactions between pair of DNA fragments that are not mapped to the genome should be included either as zeros or NAs.  
Default: FALSE.

`limit_size`  
A non-negative integer. The desired length of a log2 mean ratio vector extracted column-wise that is used at local maxima validation stage. Recommended to specify when `data_table` contains large scale genomic data (large region, chromosome or genome-wide) as it removes Hi-C bins that automatically will be removed at further stages. Ignored if NA.  
Default: NA.

`window_size`  
A non-negative integer. The desired length of Moving Average window that is used at log2 mean ratio computation stage. Recommended to specify when `data_table` contains large scale genomic data (large region, chromosome or genome-wide) as it removes Hi-C bins that automatically will be removed at further stages. Has to be specified when `limit_size` is specified. Ignored if NA.  
Default: NA.

`cores`  
A non-negative integer. The number of cores to use to run processes in parallel. If NA, apply non-parallel computing.  
Default: NA.

### Output

A list where each element named as row locus and contains a data.frame with three columns: `col_locus` is the column locus name, `col_ind` is the column locus index within the Hi-C matrix and `contact` is the interaction frequency between column locus and row locus. If the specific data.frame contains only single row (NA, NA, NA) it means that this Hi-C row contains the interactions that will be removed at further stages as the number of column positions within the row is not enough to compute the log2 mean ratios. Note that the locus names are stored in the following format `chr_start_end`. 

