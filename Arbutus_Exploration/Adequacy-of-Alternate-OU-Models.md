Adequacy of Alternate OU Models
================
Fiel Dimayacyac
02/11/2021

## Introduction

The Ornstein-Uhlenbeck model of character evolution postulates that
characters will evolve in a random walk but will be constrained by some
evolutionary optima, alpha. Brownian Motion can then be considered a
special case of OU, where the optima is 0. Alternate versions of OU
processes have been discussed in the literature, including models with
multiple evolutionary rates (OUMV), multiple adaptive optima (OUMA), and
a mixture of the two (OUMVA). By simulating OU datasets and measuring
the adequacy of each of these models to that data via the arbutus
package, I hope to gain an understanding of each of the test statistics
measured in arbutus. Measures of inadequacy for each of the alternative
models can tell me what each statistic tells me.

## Method

Simulation of phylogenetic character data was carried out via the OUwie
package. This package specializes in simulation and measurement of OU
models and variations of that model. Species trees were randomly
generated using a pure birth model via a function from the geiger
package, a package that improves and adds upon the ape package; both of
these tools are used for the analysis and handling of phylogenetic data.

Five models were tested; single rate BM, single rate and single optima
OU (labeled as OU), multirate OU (OUMV), multioptima OU (OUMA), and
multirate multioptima OU (OUMVA). 1000 sets of character data were
simulated for each model and measured for model adequacy via arbutus.
These calculations were carried out in Arbutus\_Exploration\_v2.R
script.

## Results

The main figure generated from analysis and simulation is shown below.
The figure shows the distribution and probability density of p-values
for each of the 6 test statistics generated in arbutus. High p-values
suggest that the data could be represented by the model, so a high
p-value suggests better adequacy. Each of the 6 test statistics
represents different potential violations of the data to a specific
model. Each of the models listed are the generating models, which were
then fit to a regular OU process.

![](Arbutus_Exploration/violin_all_models.png "P value distribution of test statistics")

### C.var

c.var, or the coefficient of variation, is calculated from the absolute
value of the contrasts. It accounts for rate heterogeneity. In the
figure above, it is lowest in OUMV and OUMVA, which makes sense, because
these two models have multiple rates, and therefore a regular OU model
is not accurate. Notably, it is also very low for EB, which makes sense
because the evolutionary rate is slowing down over time; so the rate is
constantly non-homogenous.

### d.cdf

This is the D statistic, which measures deviations from normality.
Because OU data is essentially random walk, it makes sense that all
models have approximately equal p-values for this statistic.
Interestingly, EB seems to deviate from normality. It would be
interesting to examine why this is.

### m.sig

This is the mean of the squared contrasts and measures if overall
evolutionary rate is over or underestimated. Because each is essentially
an OU process, violations for this statistic are not expected and not
observed. It would be interesting to see this statistic when the
parameter estimates deviate from the actual parameter values. Another
interesting area would be to look into why they are all clustered
tightly around high values. One would expect at least 5% of the p-values
to be below 0.05.

### s.asr

This calculations takes ancestral states into account to measure rate
variation in comparison to trait value. Because there are multiple rates
in OUMV and OUMVA, it follows that these two datasets do not fit the
single rate OU model. OUMA and regular OU do fit the model.

### s.hgt

This statistic measures trait variation relative to time. It seems like
all three alternative OU models are not adequate in this statistic, with
MVA being the best. Early Burst shows significant deviation for this
statistic, which makes sense, considering the slowing down of
evolutionary rate over time in the model.

### s.var

This statistic measures if contrasts deviate from expected based on
branch lengths. So if more/less evolution occurred on a shorter tree
than expected, we would see a violation. It seems that only the OU data
fit the OU model well from this measurement.

## EB Simulations

Simulating Early Burst models requires the rescaling of a tree, followed
by fitting to a Brownian Motion model. The rescaling should account for
the “slowing down” of the evolutionary rate. I first simulated a
pure-birth phylogenetic tree of 128 species, followed by rescaling the
tree with different alpha parameters. The alpha parameter represents the
strength of the EB process, i.e., how much the evolutionary rate slows
down. The original tree and the scaled trees are shown below.

![](Arbutus_Exploration/EB_tests/treeplot.png "Phylogenetic trees") Top
left represents the original tree, while top right, bottom left, and
bottom right are an alpha of half, one, and two respectively.

After using negative values for alpha, tree rescaling seems to be
working correctly as increasing alpha causes the tree to become
increasingly lengthened.

I then used fitContinuous() to test if correct alpha parameters were
being estimated. Results of the fit are as shown below.

![](Arbutus_Exploration/EB_tests/EB_fit_parameters.png) As shown in the
plot above, alpha values cluster around the values used to generate the
trees. This means the rescaling is working correctly.

## OU model parameter tests

I would also like to see what parameters in an OU model affect test
statistics generated by arbutus. I will vary alpha, sigma, and theta, as
well as test multi-regime OU models.
