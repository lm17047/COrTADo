## transform_hictable2list()

### Description

Transforms the input Hi-C table into a list where each element represents the interactions happen within single Hi-C row. The Hi-C table contains 7 columns where first three represents the interacting row locus (chromosome, start and end positions), next three represents the interacting column locus (chromosome, start and end positions) and the last one represents the interaction frequency (discrete if Hi-C data is non-normalised and continuous if Hi-C data is normalised). The table omits non-interacting pairs, i.e. if either reads of some interactions are not mapped to the genome or reads of non-interacting events, the loci pair is not represented in the table.

In case when we call TADs on a large scale genomic region (chromosome or genome-wise), we do not need the whole row to be extracted as it requires more memory and more processing time - we need only the part of the matrix that is close to the main diagonal and the size of this part is restricted by the maximum TAD width that we expect to detec. The larger is maximum TAD width, the larger are the regions at the start and the end of chromatin where TADs cannot be detected. When we work with relatively short part of chromosome, for example, it can become critical, so we want to work with whole unrestricted lower triangle matrix instead. 

### Usage

```{r}
transform_hictable2list(data_table, replace_zero = FALSE, limit_size = NA, window_size = NA, cores = NA)
```

### Required arguments

`data_table`  
A data.frame object containing 7 columns where first three represents the interacting row locus (chromosome, start and end positions), next three represents the interacting column locus (chromosome, start and end positions) and the last one represents the interaction frequency (discrete if Hi-C data is non-normalised and continuous if Hi-C data is normalised).

### Optional arguments

`replace_zero`  
Logical value indicating whether the interactions between pair of DNA fragments that are not mapped to the genome should be included either as zeros or NAs.  
Default: FALSE

`limit_size`  
A non-negative integer. The desired length of a log2 mean ratio vector extracted column-wise the is used at local maxima validation stage. Recommended to specify when `data_table` contains large scale genomic data (large region, chromosome or genome-wide) as it removes Hi-C rows that are excluded from the further analysis. Ignored if NA.  
Default: FALSE

`window_size`  
A non-negative integer. The desired length of Moving Average window that is used at log2 mean ratio computation stage. Recommended to specify when `data_table` contains large scale genomic data (large region, chromosome or genome-wide) as it removes Hi-C rows that are excluded from the further analysis. Has to be specified when `limit_size` is specified. Ignored if NA.
Default: FALSE

`cores`  
A non-negative integer. The number of cores to use to run processes in parallel. If NA, apply non-parallel computing.  
Default: FALSE
