# MethMod

## Description
MethMod can be used to perform a Gaussian Mixture Modelling of a CpG
Observed/Expected value distribution to infer the presence of DNA
methylation.

## Dependencies
MethMod was succesfully tested with these software version:

* python-3.9.2\\
* matplotlib-3.3.4
* biopython-1.78
* scikit-learn-0.24.1

Older version are likely to work as well but have not been tested.

An easy way to install the necessary tools is creating a conda
environment using the provided file 'requirements.txt':

`conda create --name <env> --file requirements.txt`

For general information how to install bioconda please visit their
website:
https://bioconda.github.io/user/install.html

## Usage
#Getting started:
`python MethMod.py --infile input.fasta`

#Getting help
`python MethMod.py -h`
usage: MethMod.py [-h] [-f] [-s] [-hd] [-w] [-sp SPECIES] -i INFILE

Gaussian Mixture modeling (GMM) to identify possible subpopulations in observed/expected (O/E) CpG distributions. The mean O/E CpG ratio for each fasta entry will be calculated. A GMM with one
or two components is performed on the distribution. A number of relevant parameters of the model are printed. Files with the distribution plots and the calculated O/E CpG ratios can be created
by using the respective option.

optional arguments:
  -h, --help            show this help message and exit
  -f, --figures         Create files with the distribution plots.
  -s, --shuffled        Perform the analysis with the input data shuffled as well.
  -hd, --header         Print a header line explaining the output columns.
  -w, --cpg             Create files containing the O/E CpG values per fasta entry.
  -sp SPECIES, --species SPECIES
                        Give a species name for you input file.
  -i INFILE, --infile INFILE
                        Input file name (fasta format).


##Example data:
cd example/;

#Simple command
`python MethMod.py --infile example/Limnephilus_lunatus.fa` 
example/Limnephilus_lunatus 	real	 0.7143 	 0.9798 	 0.2655 	 0.77 	 0.23 	 -3676.0 	 -4246.0 	 2 	 0.8 	 0.21 	 0.5

#All functions
`python MethMod.py --infile example/Limnephilus_lunatus.fa  --shuffled --figures --cpg --header`
#species	data_type	mean low	mean High	distance of means	fraction low	fraction high	AIC-1	AIC-2	Best-n	mean	stdev	skew
example/Limnephilus_lunatus 	real	 0.7124 	 0.9717 	 0.2593 	 0.76 	 0.24 	 -3676.0 	 -4244.0 	 2 	 0.8 	 0.21 	 0.5
example/Limnephilus_lunatus 	shuffled	 0.9905 	 1.0129 	 0.0223 	 0.78 	 0.22 	 -13949.0 	 -14906.0 	 2 	 1.0 	 0.14 	 0.18

The O/E CpG distribution plots of normal and shuffled data, and the
files with the O/E CpG values per fasta entry will be created in the
directory where the input file is located, i.e. example/
You can compare the files with the ones in example/exampleOutput/

NOTE: Please note, due to the nature of Gaussian Mixture Modelling the
values are varying slightly in each run. Therefore, your numbers
should be similar (first two decimal places) but not identical.

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
scikit-learn version 0.23.1 was used to estimate the parameters. The GMM was modeled with one or two
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
