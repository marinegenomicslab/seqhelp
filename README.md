
# {seqhelp}

Tidy-oriented utilities for DNA string manipulation in R

## Example usage

Fetch a subset of a fasta file (from a reference genome, for example)

``` r
library(seqhelp)

fa <- system.file("testdata", "tiny.fasta", package="seqhelp")
  
tinyseq <- fetch_seq(fasta = fa, seq_name = "seq1", seq_start = 1, seq_end = 3)

tinyseq
```

    ## [1] "ATC"

Import, manipulate, and write fasta files

``` r
tiny_fasta <- read_fasta(fa)

tiny_fasta
```

    ## # A tibble: 1 × 2
    ##   seq_id seq                           
    ##   <chr>  <chr>                         
    ## 1 seq1   ATCGATAGCGCGCGATAGCGCATGCTAGCT

``` r
tiny_fasta %>%
  mutate(new_seq = substr(seq, 1, 5)) %>%
  select(seq_id, seq = new_seq)
```

    ## # A tibble: 1 × 2
    ##   seq_id seq  
    ##   <chr>  <chr>
    ## 1 seq1   ATCGA
