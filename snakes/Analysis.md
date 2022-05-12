
# Analysis of Snake Venom Data Adequacy

------------------------------------------------------------------------

## Introduction

The data set for this analysis was taken from a
[paper](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2020.0613)
studying the expression of snake venom toxins in 52 venomous snakes. The
paper aimed to understand the varying evolutionary rates to understand
the “tempo” of evolution in adaptive radiations. They found that all
toxins undergo rate shifts, and that the Levy process better fits the
studied dynamics than either Early Burst, Ornstein Uhlenbeck, or
Brownian Motion. They concluded that there was little evidence of snake
venom exhibiting features typical of adaptive radiation, and therefore
it likely is not heavily influencing such a process. I will be analyzing
the data to test the adequacy of EB, OU, or BM models for this data.

## Summary Analysis

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.8
    ## ✓ tidyr   1.2.0     ✓ stringr 1.4.0
    ## ✓ readr   2.1.2     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
pvals <- readRDS("arbutus/pvals")

adequacy <- pvals %>% select(!m.sig) %>% transmute(c.less = c.var <= 0.05, sv.less = s.var <= 0.05, sa.less = s.asr <= 0.05, sh.less = s.hgt <= 0.05 & !is.na(s.hgt), d.less = d.cdf <= 0.05) %>% transmute(inade = c.less + sv.less + sa.less + sh.less + d.less) %>% count(inade) %>% mutate(prop = n/sum(n))

figure1 <- adequacy %>% ggplot(aes(x = inade, y = n, fill = inade)) + geom_bar(stat = "identity") + geom_text(aes(label = round(prop, digits = 2))) +
  xlab("Number of inadequacies") + ylab("Number of genes") + ggtitle("Amount of toxins by number of inadequacies in snake venom") + theme_bw() 
```

``` r
figure1
```

![](Analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

**Figure 1. Most genes show inadequacies in 3 or less test statistics.**
No genes have no inadequacies.

Overall the best-fit model is not adequate for the data.

## Results

<img src="arbutus/AIC.png" width="319" />

<img src="arbutus/arbutus.png" width="304" />

**Figure 2. Relative fit (left) and absolute fit (right) of the toxin
genes.** Most test statistics show a left-skew; suggesting inadequacies.
Most notable are c.var and s.asr; test statistics often violated when
the true model has multiple rates.

The fact that rate variation seems to be violated here when using
single-rate models, suggests that the authors of the study were correct
in using models that allow rate to change, such as the Levy “pulsed
rate” model.

Next, I will determine if a multi-rate BM model can suffice.
