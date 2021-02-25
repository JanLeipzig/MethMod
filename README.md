# MethMod

## Description
MethMod can be used to perform a Gaussian Mixture Modelling of a CpG
Observed/Expected value distribution to infer the presence of DNA
methylation.

## Dependencies
TODO

## Usage

## Parameters

## How to cite

As soon it is published please cite:
Engelhardt, et al. "Evolution of DNA methylation across Ecdysozoa."
XXX (2021): XXX

## More
Method description from Engelhardt et al. 2021:

We used two different datasets: the actual CDS data and shuffled CDS
data. For the shuffled CDS data we performed a mononucleotide
shuffling of the CDS data of each species using a custom python
script. The following analysis were performed for both the actual and
the shuffled data.

For each CDS the Observed-Expected CpG ratio was calculated using the
formula: \begin{equation} O/E_{CpG} = \frac{CG \times l}{C \times G}
\end{equation} with $C, G$, and $CG$ being the number of the
respective mono- an dinucleotids in the given CDS and $l$ being the
length of the CDS. CDS shorter than 100 nucleotides or with more than
5\% of N's in the sequence were excluded.

We used a Gaussian Mixture Model (GMM) to identify possible
subpopulations in the O/E CpG distribution. The Expectation
Maximization algorithm in the python module 'sklearn' from the library
\textit{scikit-learn} \cite{scikit-learn} version 0.23.1 was used to
estimate the parameters. The GMM was modeled with one or two
components. For the GMM with one component, we calculated the Akaike
information criterion (AIC). For the GMM with two components, we
calculated the AIC and in addition the mean of each component, the
distance $d$ of the component means and the relative amount of data
points in each component.  For the distribution of O/E CpG values, the
distribution mean, the sample standard deviation, and the skewness
were calculated as well. All pairs of parameters were analyzed using
two-dimensional scatterplots generated with R.

We used the distance between the component means as an indicator for
DNA methylation. If the distance is greater or equal to $0.25$, we
assume DNA methylation is present, otherwise it is absent.