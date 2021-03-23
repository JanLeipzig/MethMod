import numpy as np
import math
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
from sklearn.mixture import GaussianMixture
import argparse
import statistics
from scipy.stats import skew
from scipy.stats import skewtest

#Command line options
parser = argparse.ArgumentParser(description='Gaussian Mixture modeling (GMM) to identify possible subpopulations in observed/expected (O/E) CpG distributions. The mean O/E CpG ratio for each fasta entry will be calculated. A GMM with one or two components is performed on the distribution. A number of relevant parameters of the model are printed. Files with the distribution plots and the calculated O/E CpG ratios can be created by using the respective option.')

parser.add_argument("-f", "--figures", help="Create files with the distribution plots.", action="store_true")
parser.add_argument("-s", "--shuffled", help="Perform the analysis with the input data shuffled as well.",action="store_true")
parser.add_argument("-hd", "--header", help="Print a header line explaining the output columns.", action="store_true")
parser.add_argument("-w", "--cpg", help="Create files containing the O/E CpG values per fasta entry.",action="store_true")
parser.add_argument("-sp", "--species", help="Give a species name for you input file.")
parser.add_argument("-i", "--infile", help="Input file name (fasta format).", required=True)

args = parser.parse_args()

infile=args.infile
if args.species:
    species=args.species
else:
    species=infile.split(".")[0]

if args.header:
    print("#species\tdata_type\tmean low\tmean High\tdistance of means\tfraction low\tfraction high\tAIC-1\tAIC-2\tBest-n\tmean\tstdev\tskew")


#Custom methods
def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier

def get_best_n(lines_local):
    X_local=lines_local.reshape(-1,1)

    gmm1_local = GaussianMixture(n_components=1)
    label1_local=gmm1_local.fit(X_local).predict(X_local)
    aic1_local=gmm1_local.aic(X_local)

    gmm2_local = GaussianMixture(n_components=2)
    label2_local=gmm2_local.fit(X_local).predict(X_local)
    aic2_local=gmm2_local.aic(X_local)

    best_n_local=2
    if aic1_local < aic2_local:
        best_n_local=1
    return best_n_local


##Normal CpG rates
#Calculate O/E CpG rate start

mylist = []

for record in SeqIO.parse(infile, "fasta"):
    g=record.seq.upper().count("G")
    c=record.seq.upper().count("C")
    cg=record.seq.upper().count("CG")
    n=record.seq.upper().count("N")

    l=len(record.seq)

    if (c*g)>0:
        oe=(cg*l)/(c*g)
    
        if l>=100 and (n/l)<=0.05 and oe>0:
            #print(oe)
           mylist.append(oe)

if args.cpg:
    with open(species+'.cpg', 'w') as f:
        for item in mylist:
            f.write("%s\n" % item)
#Calculate O/E CpG rate end


#Gaussian Mixture Modelling start
lines=np.array(mylist)
X=lines.reshape(-1,1)

gmm1 = GaussianMixture(n_components=1)
label1=gmm1.fit(X).predict(X)
aic1=gmm1.aic(X)

gmm2 = GaussianMixture(n_components=2)
label2=gmm2.fit(X).predict(X)
aic2=gmm2.aic(X)

best_n=get_best_n(lines)

if float(gmm2.means_[0])<float(gmm2.means_[1]):
    meanLow=0
    meanHigh=1
else:
    meanLow=1
    meanHigh=0

lowPoints=lines[label2==meanLow]
highPoints=lines[label2==meanHigh]
fractionLow=round_half_up(len(lowPoints)/len(lines),2)
fractionHigh=round_half_up(len(highPoints)/len(lines),2)

meanLowRound=round_half_up(gmm2.means_[meanLow],4)
meanHighRound=round_half_up(gmm2.means_[meanHigh],4)
meanDist=round_half_up(gmm2.means_[meanHigh]-gmm2.means_[meanLow],4)

stdev=round_half_up(statistics.stdev(lines),2)
mean=round_half_up(statistics.mean(lines),2)
skew_value=round_half_up(skew(lines),2)

#print(species,"real\t",meanLowRound,"\t",meanHighRound,"\t",meanDist,"\t",fractionLow,"\t",fractionHigh,"\t",round_half_up(aic1,0),"\t",round_half_up(aic2,0),"\t",best_n)
print(species,"\treal\t",meanLowRound,"\t",meanHighRound,"\t",meanDist,"\t",fractionLow,"\t",fractionHigh,"\t",round_half_up(aic1,0),"\t",round_half_up(aic2,0),"\t",best_n,"\t",mean,"\t",stdev,"\t",skew_value)

#Gaussian Mixture Modelling end



##Normal Data
if args.figures:
    plt_lin=lines
    M_best = gmm2

    bins=100
    x = np.linspace(0, 2, 1000)
    logprob = M_best.score_samples(x.reshape(-1, 1))
    responsibilities = M_best.predict_proba(x.reshape(-1, 1))
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]
    
    plt.ylim(top = 4, bottom = 0)
    plt.hist(plt_lin, bins, density=True, histtype='stepfilled', alpha=0.4) #original data
    plt.plot(x, pdf, '-k')
    plt.plot(x, pdf_individual, '--k')
    plt.ylabel('Density')
    plt.xlabel('O/E CpG')
    plt.title(species)
    
    plt.savefig(species+"-fit.png",dpi=800)




##Shuffled CpG rates
if args.shuffled:

    #Calculate O/E CpG rate start
    mylist_shuf = []
    
    for record in SeqIO.parse(infile, "fasta"):

        nuc_list=list(record.seq)
        random.shuffle(nuc_list)
        shuffled_rec=SeqRecord("".join(nuc_list),'','','')
        g=shuffled_rec.seq.upper().count("G")
        c=shuffled_rec.seq.upper().count("C")
        cg=shuffled_rec.seq.upper().count("CG")
        n=shuffled_rec.seq.upper().count("N")        
        l=len(shuffled_rec.seq)

        if (c*g)>0:
            oe=(cg*l)/(c*g)
            
            if l>=100 and (n/l)<=0.05 and oe>0:
                mylist_shuf.append(oe)
    if args.cpg:
        with open(species+'.shuffled.cpg', 'w') as f:
            for item in mylist_shuf:
                f.write("%s\n" % item)
    #Calculate O/E CpG rate start

    #Gaussian Mixture Modelling start
    lines_shuf=np.array(mylist_shuf)
    X_shuf=lines_shuf.reshape(-1,1)
    best_n_shuf=get_best_n(lines_shuf)

    gmm1_shuf = GaussianMixture(n_components=1)
    label1_shuf=gmm1_shuf.fit(X_shuf).predict(X_shuf)
    aic1_shuf=gmm1_shuf.aic(X_shuf)

    gmm2_shuf = GaussianMixture(n_components=2)
    label2_shuf=gmm2_shuf.fit(X_shuf).predict(X_shuf)
    aic2_shuf=gmm2_shuf.aic(X_shuf)

    if float(gmm2_shuf.means_[0])<float(gmm2_shuf.means_[1]):
        meanLow_shuf=0
        meanHigh_shuf=1
    else:
        meanLow_shuf=1
        meanHigh_shuf=0

    lowPoints_shuf=lines_shuf[label2_shuf==meanLow_shuf]
    highPoints_shuf=lines_shuf[label2_shuf==meanHigh_shuf]
    fractionLow_shuf=round_half_up(len(lowPoints_shuf)/len(lines_shuf),2)
    fractionHigh_shuf=round_half_up(len(highPoints_shuf)/len(lines_shuf),2)
    
    meanLowRound_shuf=round_half_up(gmm2_shuf.means_[meanLow_shuf],4)
    meanHighRound_shuf=round_half_up(gmm2_shuf.means_[meanHigh_shuf],4)
    meanDist_shuf=round_half_up(gmm2_shuf.means_[meanHigh_shuf]-gmm2_shuf.means_[meanLow_shuf],4)

    stdev_shuf=round_half_up(statistics.stdev(lines_shuf),2)
    mean_shuf=round_half_up(statistics.mean(lines_shuf),2)
    #skew_shuf=round_half_up(skew(lines_shuf),2)
    skew_value_shuf=round_half_up(skew(lines_shuf),2)

    print(species,"\tshuffled\t",meanLowRound_shuf,"\t",meanHighRound_shuf,"\t",meanDist_shuf,"\t",fractionLow_shuf,"\t",fractionHigh_shuf,"\t",round_half_up(aic1_shuf,0),"\t",round_half_up(aic2_shuf,0),"\t",best_n_shuf,"\t",mean_shuf,"\t",stdev_shuf,"\t",skew_value_shuf)

    
    #print(species,"shuffled\t",meanLowRound_shuf,"\t",meanHighRound_shuf,"\t",meanDist_shuf,"\t",fractionLow_shuf,"\t",fractionHigh_shuf,"\t",round_half_up(aic1_shuf,0),"\t",round_half_up(aic2_shuf,0),"\t",best_n_shuf)
    #Gaussian Mixture Modelling end

    #Plot histograms start
    if args.figures:
        plt.clf()
        plt_lin=lines_shuf
        M_best_shuf = gmm2_shuf
            
        logprob = M_best_shuf.score_samples(x.reshape(-1, 1))
        responsibilities = M_best_shuf.predict_proba(x.reshape(-1, 1))
        pdf = np.exp(logprob)
        pdf_individual = responsibilities * pdf[:, np.newaxis]
        
        plt.ylim(top = 4, bottom = 0)
        plt.hist(plt_lin, bins, density=True, histtype='stepfilled', alpha=0.4) #original data
        plt.plot(x, pdf, '-k')
        plt.plot(x, pdf_individual, '--k')
        plt.ylabel('Density')
        plt.xlabel('O/E CpG')
        plt.title(species)
            
        plt.savefig(species+"-fit.shuffled.png",dpi=800)
    #Plot histograms end
