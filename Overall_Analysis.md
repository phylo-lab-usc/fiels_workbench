
# Analysis of the Adequacy of Widely Used Phylogenetic Models in Gene Expression Data

------------------------------------------------------------------------

## Introduction

The analyses in this document were run on 9 different gene expression
data sets, outlined in the README file in this repository. These data
sets are of hundreds of different species and tens of thousands of
different genes. Some analyses were run using the species-level
phylogeny, while others were run using individual gene-level
phylogenies. The models interrogated in this analysis include the
Brownian Motion model depicting random walk of a character trait, the
Ornstein-Uhlenbeck model depicting random walk with a constraint value,
and the Early Burst model depicting random walk with a slowing down of
evolutionary rate as time passes. A fourth model is added, Brownian
Motion with multiple rates, that allows different branches of the tree
to have separate evolutionary rate parameters but is otherwise Brownian
Motion. In the literature, each of these models is often used as a
stand-in of an evolutionary process; Brownian Motion as neutral drift,
Ornstein-Uhlenbeck as stabilizing selection, and Early Burst as adaptive
radiation. The selection of a model as the “dominant model” through
various model-selection processes is then often used to state that
stabilizing selection is the dominant mode of evolution for some set of
genes for example. Rather than make specific claims about the results of
each paper, this analysis aims to understand and highlight patterns of
model adequacy in the field; i.e., do the models we are using actually
describe the variation in the data? Or is there much left to be desired?

## Summary

By using the Akaike Information Criterion to measure relative fit, I
show that broadly, the Ornstein-Uhlenbeck model is the best fit for most
of the data; i.e., the model with the lowest AIC for each gene is most
often OU by a fair margin. In the absolute sense, most of the data was
adequately modeled by the best-fit model (most commonly OU), but these
models were very often inadequate in terms of rate heterogeneity. I then
ran simulations to show that processes that truly had multiple rate
parameters would violate the same test statistics as seen in the data.
To solve this I then added a fourth model to the analysis; multi-rate
Brownian Motion. In a relative sense, addition of the model then shifted
the domination to this multirate BM so that most of the data was
best-fit by a multi-rate model. In an absolute sense, it did somewhat
reduce issues with rate heterogeneity, but other issues still persist.
All in all, it seems that in general, rate heterogeneity is not being
accounted for in the most commonly used models of evolution in the
field, and usage of those types of models would improve data capture.

## Data used in this analysis

| Data set                      | Number of Genes | Number of Species |
|-------------------------------|-----------------|-------------------|
| cichlids                      | 32,596          | 73                |
| fishes                        | 16,740          | 20                |
| Heliconius\_Butterflies       | 2,393           | 5                 |
| Mammal\_organs                | 5,320           | 10                |
| snakes                        | 11              | 52                |
| comparative\_expression\_2017 | 8,333           | 9                 |
| GeneExpression\_coevolution   | 3,556           | 18                |
| interspecific\_rnaseq         | 3,560           | 14                |
| amalgam data                  | 1,377           | 21                |

## Simulations Illustrate Expected Inadequacy Patterns

To visualize patterns I expect to see when data is adequate or
inadequate, I first performed a set of simulations. In these
simulations, the generating model varies, but the model being fit to the
data stays consistent. In this way, I can determine where I can expect a
model to be inadequate when the true model does not align with the model
we are trying to fit. Specifically, the model being fit is an
Ornstein-Uhlenbeck process; the model most commonly used in Phylogenetic
Comparative Methods. The results of these simulations are depicted
below.

![](Arbutus_Exploration/violin_all_models.png)

**Figure 1. Model Adequacy when a series of processes are fit to a model
that is not the true one.** P values are depicted on the y-axis and each
test statistic is labeled. The generating models are listed in the
legend on the right, where models starting with “M” are variants of an
OU process with multiple optima or multiple alpha values.

The models tested in these simulations include BM, BMS (multiple rate
brownian motion), EB, OU, and multiple optima OU variants with either
multiple alpha values (MA), evolutionary rates (MV), or both (MVA).

C.var is the coefficient of variation, used to test if rate
heterogeneity is accounted for by a fitted model. As shown by the
simulation, a process that is truly Early Burst will show stark
inadequacies in this test statistic; notice the high build-up of
p-values near 0.05.

## Support For Currently Used Models

Based on my analysis there is evidence supporting the usage of the
models mentioned above, specifically the OU model. The amalgamated data
set is the largest data set showing support for the models above. This
is an amalgamation of 1903 RNA Seq studies with 6 organs and 21
vertebrate species. In a relative sense, the data was best fit by the
Ornstein-Uhlenbeck model and in the absolute sense, test statistics show
wide distributions; suggesting adequacy of the model in all 5 aspects.
The adequacy analyses are shown below.

![A](Overall_Analysus_Plots/figure2a.png)

![](Overall_Analysus_Plots/figure2b.png)

**Figure 2. Relative Fit of the three tested models (top) and Absolute
Fit of the best-fit model (bottom) for Amalgamated Data.** The bulk of
the data was best fit by an Ornstein-Uhlenbeck process. The adequacy of
the best fit model for each gene was high across the tested statistics,
with evenly distributed test statistics across the board.

As shown by Figure 2, the amalgamated data set showed high adequacy for
the best-fit model, which tended to be the Ornstein-Uhlenbeck model.
This supports the current paradigm of PCMs across different species and
thousands of genes. In fact, 77% of genes were shown to be fully
adequate across all statistics, which increases to 96% when including
genes with only one inadequacy.

One of the goals of this analysis was to identify models that may
describe gene expression data better in the context of phylogenetic
comparative analysis. In the “Coevolution” data set, the author used a
BM model to identify evidence of coevolution of proteins across 18 fungi
species. I then performed arbutus analysis of the expression data after
fitting to a BM model to see how well the model explains the variation
in the data. A summary of the analysis is shown below.

![](Overall_Analysus_Plots/figure3.png)

**Figure 3. Number of genes by number of inadequacies (cut off of 0.05)
across test statistics for Coevolution data set.** Just over half of the
data is adequately modeled by a BM process, with some of the genes being
inadequate in all test statistics.

I then performed relative fit analysis to see if there were models that
explained the data better, with my results shown below.

![](Overall_Analysus_Plots/figure4.png)

**Figure 4. Count of the model with the lowest AIC value for each gene
in the Coevolution data set.** The OU model is overwhelmingly supported
by the data.

Finally, I performed arbutus analysis on the best fit model for each
gene.

![](Overall_Analysus_Plots/figure5.png)

**Figure 5. Number of genes by number of inadequacies (cut off of 0.05)
across test statistics for Coevolution data set when using the best fit
model.** Genes across the board show great decreases in overall
inadequacy.

## Species Phylogenies Are More Adequate Than Gene Family Phylogenies

One interesting finding from the analysis of the Coevolution data set
was just how adequate the data was compared to data analyzed in previous
studies. A previous analysis found current models to be quite inadequate
to describe gene expression data. This analysis is recreated in the
Comparative Expression repository. One major difference between these
studies was the usage of a species phylogenetic tree in the case of the
Coevolution data set, in contrast to gene-family phylogenies in the
comparative expression data. To uncover if this adequacy difference was
due to the difference in methodologies; I, along with Doris Wu,
re-analyzed the Coevolution data set substituting the species
relationships with gene-family phylogenies, so that each group of genes
was described by the relationship between those genes. The generation of
the gene family phylogenetic trees is described in this
[repository](https://github.com/pennell-lab-ubc/gene-phylogeny-pipeline).
The summary of the arbutus analysis of this data set is shown below.

![](Overall_Analysus_Plots/figure6.png)

**Figure 6. Inadequacies Found When Using Gene Family Phylogenies.**
There are more inadequate genes across the board, but less NA values.

This comparison was also carried out in the Mammal Organs data set. The
p-value distributions for the test statistics are displayed below.

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 4297 rows containing non-finite values (stat_bin).

![](Overall_Analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

**Figure 7. Comparison of Test Statistic Distributions When Using Gene
Family Phylogenies vs Species Phylogenies in Analyzing Adequacy of
Mammal Organs Expression Data.** The inadequacies seen in c.var and
s.asr appear to be higher when using species phylogenetic relationships,
but inadequacies in s.var and s.hgt increase when using gene family
phylogenies.

As shown by the figure above, inadequacies increase in certain aspects,
but decrease in others. Interestingly, using a species phylogeny seems
to cause many genes to show no phylogenetic signal; shown by the drop
off in p values near the low end for s.hgt, which was also seen in the
Coevolution data set. This was confirmed in the Coevolution data set to
be due to low phylogenetic signal in the plot below.

![](Overall_Analysus_Plots/figure8.png)

**Figure 8. Phylogenetic Signal of Genes with NA Values in the S.hgt
statistic have Low Phylogenetic Signal.** The K value describes
phylogenetic signal where higher values indicate higher signal and the
p-values are the result of a hypothesis test where p-values of 0.05 are
genes with phylogenetic signal.

In contrast, using gene family phylogenies in analysis causes genes that
have no phylogenetic signal when using a species phylogeny to have
signal, but they are inadequate in the S.hgt statistic. This statistic
shows inadequacies when trait variation over time is not accounted for.

Together, it seems that using gene family phylogenies over species
phylogenies increases the amount of inadequacies in a set of data, but
this increase may be only due to NA values showing up as inadequate
p-values. This suggests that using a species phylogeny can allow the
best-fit model to describe the data better statistically (for genes not
thrown out for NA p values), but at the expense of the phylogenetic
signal of many genes, which may decrease the ability to make biological
inferences. In essence, because the species phylogeny is an average and
summation of the relationships between all the genes, it is often a
“close enough” tree, but for genes that do not fully follow the
relationship at the species level, the model that describes their
evolution may be misleading.

## Many Data Sets Show Model Inadequacy For Rate Heterogeneity

While the best-fit model for some data sets has been shown to be highly
adequate, many of the other data sets in this analysis showed similar
patterns in where the best-fit models tended to fail.

![](Overall_Analysus_Plots/figure9a.png)

![](Overall_Analysus_Plots/figure9b.png)

**Figure 9. P-value Distributions of Arbutus Test Statistics for
Extremophile Fishes Data Set (top) and Cichlids Data Set (bottom).**
Inadequacies are seen in all test statistics except for d.cdf, but most
notably in c.var and s.asr for extremophile fishes, whereas the cichlids
data set shows inadequacies in mainly c.var, d.cdf, and s.asr.

In general, the most commonly violated test statistics are the
coefficient of variation (c.var) and the statistic of ancestral state
reconstruction (s.asr) as shown by Figure 7 and 9. These test statistics
both tend to be violated together, and are both in some way related to
heterogeneity in evolutionary rate over the tree. Going back to Figure
1, these test statistics are both violated when the true model has
multiple rates, but the fitted model does not (MV and MVA). The simplest
multi-rate model is the BMS, or multiple rate BM, so I ran simulations
to see what p-value distributions would look like if the true model were
multirate BM, and the other, single-rate models were fit.

![](Arbutus_Exploration/Figures/multirate_comparison.png)

**Figure 10. P value Distributions of Each Test Statistic When the True
Model is BMS.** The fitted models are listed in the legend to the side.
Both C.var and S.asr show clusters of low pvalues when the fitted model
is any other model than BMS.

As shown by the simulation, fitting models with one evolutionary rate
regime (BM, OU, EB) to a truly multi-rate process would also generate
high inadequacies in the C.var and S.asr test statistics, confirming
this as a possible cause of inadequacies in data. I then tested the
addition of a BMS process to the Arbutus Analysis process on multiple
data sets.

![](Overall_Analysus_Plots/figure11.png)

**Figure 11. Addition of the BMS model lowers adequacy violations in
C.var and S.asr in Extremophile Fish data.** More violations are seen in
the s.hgt statistic and s.var statistic when a BMS model is included in
analysis.

![](Overall_Analysus_Plots/figure12.png)

![](Overall_Analysus_Plots/figure12b.png)

**Figure 12. Addition of a BMS Model to Analysis Increases the Amount of
Genes with No Inadequacies in Extremophile Fish Data.** Genes with 3 or
more inadequacies showed little change with the addition of a BMS model.

This same pattern of increasing the proportion of genes with fully
adequate models by mainly decreasing inadequacies in C.var and S.asr was
also seen in the Interspecific RNA Seq data set.

One issue with using a model with multiple evolutionary rates is the
need to pre-define selective regimes. This means that one would have to
have some biological question or knowledge specific to a set of taxa.
However, when broadly analyzing the adequacy of models themselves over
the field, this leaves the problem of how to define evolutionary regimes
to allow for comparison between data sets. In cases where the authors of
the original study made inferences about evolutionary regimes, I simply
used those. However, for data sets such as the Snake Venom and Cichlids
data sets, I instead used the Motmot tool. This tool identifies clades
that experience a shift in evolutionary rate through multiple hypothesis
testing. The results of using motmot combined with Arbutus are shown
below.

![](Overall_Analysus_Plots/figure13.png)

![](Overall_Analysus_Plots/figure13b.png)

**Figure 13. Adequacy of the Best Fit Model for Snake Venom Data (top)
and Cichlids Data (bottom) with and without inclusion of BMS models to
Analysis.** Adequacy increased for snake venom genes, but no change was
seen in cichlids data.

Overall, the addition of multi-rate BM models to arbutus analysis did
increase the adequacy of the best-fit model for the data. Specifically,
it decreased issues with C.var and S.asr, which both detect violations
when there are issues with evolutionary rate heterogeneity across a
tree. However, the magnitude of this increase ranges from moderate
(Figure 12) to very little (Figure 13). Thus, factors other than rate
heterogeneity seem to be coming into play that is causing models to be
inadequate even with the addition of a model that allows for multiple
evolutionary rates. Early Burst models were shown to be violated in this
test statistic as well (Figure 1), suggesting that one of these factors
may be “static” evolutionary rates that are seen in all the models used
except for EB.

## Conclusion

Model adequacy is the ability of a model to explain the variation seen
in a set of data. More adequate models allow researchers to make better
predictions about some statistical process. In this case, the models we
are using are proxies for evolutionary processes such as neutral drift,
adaptive selection, divergent evolution, etc. Phylogenetic Comparative
Methods rely on such models to make claims about the evolution of
organisms, but these methods have not been extensively evaluated in the
context of gene expression data. In this analysis I have shown that of
the three models widely used in the field today, the Ornstein-Uhlenbeck
model tends to be the best-fit model for much of the data in a relative
sense. In an absolute sense, the results are mixed. For data sets that
use a species phylogeny as the basis, the model does a relatively
adequate job. While this is good news statistically, the adequacy of
this model tends to decrease a fair amount when a more accurate gene
family phylogeny is used for each gene studied. This suggests that using
a species phylogeny may allow the current models to perform well
statistically, it may decrease the predictive power or ability for
biological inference. Additionally, when these models do fail, they tend
to be related to issues with heterogeneity in evolutionary rate.
Simulation data showed that fitting data that is truly multi-rate with a
single rate model recreated these patterns. Addition of a multi-rate
model to analysis did somewhat alleviate the inadequacies with rate
heterogeneity, but it did not eliminate them. Therefore, other factors
may be at play. These factors may include differences in gene expression
normalization, static evolutionary rates, batch effects, etc. I suggest
that models that account for intraspecific variation, such as the
EvoGeneX modified OU model, may further alleviate such model
inadequacies, allowing researchers to make better inferences on
evolutionary data.
