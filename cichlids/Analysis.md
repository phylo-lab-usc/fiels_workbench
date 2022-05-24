
# Analysis of African Cichlid Fishes Data Adequacy

------------------------------------------------------------------------

## Introduction

The data set is taken from a
[paper](https://www.nature.com/articles/s41559-020-01354-3) studying
gene expression data in six different organs for 73 species of African
Cichlid found in Lake Tanganyika. Specifically, the paper analyzes the
expression patterns of both protein-coding genes and lncRNA to study the
evolutionary dynamics associated with the rapid adaptive radiation known
from these fish.They find that the rate of gene expression evolution
varies between organs, and that the noncoding transcripts (lncRNA)
evolve faster than coding transcripts. Interestingly, they also found
that the rate of evolution accelerated later rather than earlier.
Lastly, they used model fitting between EB, OU, and BM data to suggest
that most of this evolution was dominated by stabilizing selection, with
OU acting as the proxy for stabilizing selection. I will be analyzing
this data set to see if that last point stands by testing the adequacy
of each of those models for this data set.

## Summary Analysis

![](Analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

**Figure 1: Using the best-fit model (chosen by AIC) shows a very low
adequacy overall.** A) Total adequacy for long-non-coding RNA. B) Total
adequacy for protein-coding genes. Adequacy overall is very low, with
genes being only 4% adequate for lncRNA and 24% for protein-coding
genes.

Overall it seems like the best-fit models did not adequately capture the
data in neither lncRNA nor protein coding genes, suggesting that other
model(s) may be needed.

## Results

### Initial Arbutus Analysis

![](Analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 696 rows containing non-finite values (stat_bin).

![](Analysis_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

**Figure 2. Relative fit (A) and absolute fit (B) of the protein-coding
genes.** Overall, a OU model fits the data the best in a relative sense,
and in an absolute sense the best-fit model is quite inadequate. C.var,
d.cdf, and s.asr show meaningful inadequacy.

<img src="arbutus/AIC/AIC_br_lnc.png" width="327"/>

<img src="arbutus/figures/arbutus_br_lnc.png" width="326"/>

**Figure 3. Relative fit (left) and absolute fit (right) of the lncRNA
genes.** Overall, a OU model fits the data the best in a relative sense,
and in an absolute sense the best-fit model is very inadequate. C.var,
d.cdf, s.asr, and s.var show meaningful inadequacy.

Similarly to what the original paper found, most genes seemed to fit an
OU model the best for both protein-coding and lncRNA genes. However, in
an absolute sense the best fit models had high and very high numbers of
inadequacies to protein-coding and lncRNA genes respectively. Because
these inadequacies were related to c.var and s.asr, I hypothesize that a
muti-rate BM model will better describe the data.

However, to perform fitting using a BMS model, selective regimes need to
be defined. To get around that issue I will use the
[Motmot](https://github.com/PuttickMacroevolution/motmot) package to
identify rate shifts. This package works by testing likelihoods of a
rate shift per clade. Likelihood of a rate shift higher than a certain
cutoff value will be considered shifts in evolutionary rate for that
gene. The tree will then be split by evolutionary regimes before and
after the rate shift; taking the guesswork out of the equation.

![](Analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

**Figure 4. Motmot analysis finds one rate shift in the tree.**

Adequacy analysis of the protein data using a multiple-rate BM model
produces the following figures.

![](Analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 890 rows containing non-finite values (stat_bin).

![](Analysis_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

**Figure 5. Relative fit (A) and absolute fit (B) of protein-coding
genes with addition of a BMS model.** Overall, a OU model fits the data
the best in a relative sense, with BMS being a close second. In an
absolute sense the best-fit model is very inadequate. C.var, d.cdf,
s.asr, and c.var show meaningful inadequacy even while using a BMS
model.

Inadequacy of the data remains high even while using a BMS model,
showing other factors may be at play causing inadequacies in c.var. I
then wanted to quantify the effect of the addition of BMS.

![](Analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

**Figure 6: Addition of a multi-rate BM process does little to lower the
inadequacies found in the data, with similar proportions of inadequacies
before and after usage.**

Overall this shows that with this data set there must be other
considerations not due to multiple rates on different branches.
