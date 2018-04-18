# PICS-Ord: pairwise identity and cost scores ordination
Extract phylogenetic information from ambiguously aligned regions

## Description

PICS-Ord is an algorithm to extract phylogenetic information from
hard-to-align regions of multiple-sequence alignments. It has been
implemented in an R-based program using [Ngila](https://github.com/reedacartwright/ngila/) as a back end.

### Citation

Lücking R, Hodkinson B, Stamatakis A, and Cartwright RA (2011) PICS-Ord: unlimited coding of ambiguous
regions by pairwise identity and cost scores ordination. BMC Bioinformatics, 12:10. doi:[10.1186/1471-2105-12-10](https://dx.doi.org/10.1186/1471-2105-12-10).

URL: <https://github.com/CartwrightLab/PICS-Ord/>

## Getting Started

In order to run PICS-Ord, you need to install **R** on your computer.
Download the most recent version of R from <http://www.r-project.org/>
and follow the installation instructions.

You will also need to download and install **Ngila Release-1.3** or
greater from <https://github.com/reedacartwright/ngila/releases>.

After downloading the PICS-Ord package, unzip the file folder and store
it on your desktop or in your preferred directory. The folder contains
the application and documentation files.

## Terminology Convention

By "ambiguous regions", we refer to portions of a multiple fixed,
sequence alignment (MSA) that are aligned with low confidence. These are
regions where different alignment configurations yield the same
alignment score or the optimal alignment does not have a significantly
higher score than many sub-optimal alignments.

## Preparing Input Files

PICS-Ord accepts an ambiguous region of the MSA provided as a
Fasta-formatted file (suffix .fas or .fasta), Phylip-formatted file
(suffix .phy), or Clustal-formatted file (suffix .aln). It then produces
a Phylip-formatted output file. The output file contains an "alignment"
of encoded characters. Note that, you need to identify the ambiguous
regions in your MSA in your own, prior to invoking PICS-Ord, and save
each such region as separate fasta file. We recommend to first align
your sequence using MAFFT, <http://mafft.cbrc.jp/alignment/software/>.

After alignment with MAFFT, you will need to delimit/identify the
ambiguous regions manually. This can be done visually by looking for
regions with a high gap or substitution content. However, it is more
informative to use a tool devoted to this task. We recommend running the
unaligned sequence assembly through GUIDANCE at
<http://guidance.tau.ac.il/> using the default settings.

GUIDANCE will run your sequence assembly through MAFFT and assign
alignment confidence scores to each column; among other results, it will
return a html file named "MSA Colored" according to the confidence
score. In this graphical file, sites aligned with high confidence appear
in shades of magenta and pink, whereas sites aligned with low
confidence, that probably correspond to ambiguous regions, appear in
shades of blue. You may then compare the results of your initial manual
identification of ambiguous regions and adjust it as required.

Based on our experience, GUIDANCE generally produces good results in
delimiting ambiguous regions, but can contain minor alignment errors
that need to be adjusted manually; therefore our recommendation is to
use both manual/visual delimitation of ambiguous regions prior and
GUIDANCE.

As already mentioned, after having identified ambiguous regions, you
need to save each region (if there is more than one ambiguous region in
the MSA) as a separate sequence file. To do so, you need to save the
entire MSA under different names for each ambiguous region that you
intend to recode with PICS-Ord. Then, open each file and delete the
alignment columns before and after the specific ambiguous region.

*It is not necessary to include conserved flanking regions in the
ambiguous region and you should not do so if you intend to use the
free-end-gaps option in NGILA (see below).*

If you have not placed `picsord` in your path, you will need to copy
these input files to whatever directory you placed `picsord`.

## Recoding Ambiguous Regions Using PICS-Ord

PICS-Ord runs from the command line. The PICS-Ord command line is

```
$ Rscript --vanilla picsord.R input.fas > output.phy
```

where `input.fas` is your ambiguous region input file and `output.phy`
is the recoded output file. For multiple ambiguous regions, you will
need to invoke this command for each ambiguous region fasta input file,
e.g.

```
$ Rscript --vanilla picsord.R inputRegion1.fas > outputRegion1.phy
$ Rscript --vanilla picsord.R inputRegion2.fas > outputRegion2.phy
$ Rscript --vanilla picsord.R inputRegion1.fas > outputRegion1.phy
```

The application will return an individual Phylip file for each ambiguous
region. The number of columns/sites in those Phylip files corresponds to
the number of ordination axes with positive eigenvalues, minus the
number of axes that contain all-zeros after rescaling and recoding.

## PICS-Ord Program Parameters

PICS-Ord uses three modules to compute ambiguous region codes:

 1. computing a pairwise distance matrix using Ngila
 1. ordinating the distance matrix by means of principal coordinate analysis (PCoA)
 1. rescaling and recoding the ordination axis scores to integer values

Ngila program parameters can be adjusted, as well as the parameters for
rescaling and recoding the ordination axis scores. The latter is not
recommended, since parameters are optimized for producing simple-digit
integer scores between 0 and 9 that can seamlessly be analyzed within a
phylogenetic context.

The default Ngila command line parameters in the `picsord.R` script are

```
ngilacmds <- "-m zeta -o dist-c:-"
```

There is also second line that tells Ngila to compute pairwise distance
scores for all pairs of sequences. Since this is required to produce the
complete distance matrix, this command line should not be changed.

The first line contains the main Ngila parameters: the underlying model
(-m) for computing the cost scores (distances) and additional parameters
such as the output file name (-o). A detailed explanation of parameter
settings is available at
<https://github.com/reedacartwright/ngila>.

A user will probably be most interested in changing the pair-wise
distance score model and in using the free end gaps option. For
computing distances, five models are available: zeta, geo, aazeta,
aageo, and cost. Zeta and geo are models that calculate alignment costs
based on evolutionary parameters. They are the same except that the
former uses a biologically more realistic power-law distribution for
indel lengths, while the latter uses the more common geometric
distribution. (Aazeta and aageo are equivalent models for proteins.) The
cost model allows you to directly specify alignment costs. Reasonable
modifications of the `picsord.R` command line may look as indicated
below:

```
ngilacmds <- "-m zeta -o dist-c:-"
ngilacmds <- "-m geo -o dist-c:-"
ngilacmds <- "-m cost -o dist-c:-"
```

The free end gaps option \[\--free-end-gaps\] allows gaps at the start
and end of the smaller sequence to have lower or no cost than gaps
within the sequence. This is particularly useful when it is expected
that the end points of a sequence pair are not homologous. This will be
especially useful for ambiguous regions in ribosomal DNA that correspond
to transcribed, but later degraded spacer regions. We found that,
applying this option results in ambiguous region codes that are more
homogeneous within subtrees and hence increase bipartition support. The
corresponding `picsord.R` command line using, for instance, a zeta
alignment cost function would look like this:

```
ngilacmds <- "-m zeta --free-end-gaps -o dist-c:-"
```

## Integration with Phylogenetic Analysis

The recoded output files in Phylip format can easily be concatenated
with the remaining non-ambiguous MSA either manually or using MESQUITE
(<http://mesquiteproject.org/mesquite/download/download.html>). Analysis
under parsimony is straightforward, and we recommend using the PICS-Ord
codes as ordered characters. If you intend to use unordered characters,
you should weight the ambiguous region columns according to the range of
values encountered in each column; e.g. a column with the maximum value
7 and the minimum value 0 receives the weight 7, whereas a column with
the maximum value 2 and the minimum value 1 receives the weight 1. The
latter is useful when analyzing the dataset under maximum likelihood.

For analysis under maximum likelihood, we recommend using RAxML version
7.2.6 or later (the most recent version is v728 ALPHA, which is already
very stable), which is available at
<https://sco.h-its.org/exelixis/web/software/raxml/>. Please Follow the
instructions in the manual and the on-line help (by typing
`./raxmlHPC -h`) for compiling and installing RAxML. Below, we provide
a simple example for analyzing a partitioned MSA containing DNA data and
PICS-Ord recoded partitions with RAxML.

Consider a PHYLIP-formatted DNA MSA with two PICS-Ord re-coded regions
that may look as follows:

```
4 14
Seq1 0123ACGTTT2345
Seq2 1022ACGGTT2234
Seq3 0031AGTTTT1201
Seq4 1010ACGTTT2101
```

To analyze this data set we will first need to create a partition file
that tells RAxML how to partition the data. This is just a plain text
file that will then need to be passed to RAxML. This file would look as
follows in the concrete example:

```
MULTI, p1 = 1-4
DNA, p2 = 5-10
MULTI, p3=11-14
```

We may save this information, e.g. as `part.txt`. This file tells RAxML
that sites 1--4 of the alignment are multi-state sites and the partition
will be named "p1". (You can chose an arbitrary name here.) In fact,
every line of this file describes one partition of the dataset. In the
second line we tell RAxML that we have a DNA partition that stretches
from sites 5 to 10, which we name p2. And finally, the last partition
(the second ambiguous region) is once again a multi-state partition that
encompasses sites 11--14.

Given the two files, we can then do, e.g., a standard RAxML ML tree
search by typing:

```
./raxmlHPC -q part.txt -s picsOrd.phy -m GTRCAT -n T1
```

This will conduct one tree search under the CAT model of rate
heterogeneity for each partition. If we want to use GAMMA we should
type:

```
./raxmlHPC -q part.txt -s picsOrd.phy -m GTRGAMMA -n T2
```

The two above commands will use a GTR substitution model for the two
multi-state (ambiguous) partitions by default. If you want to use, e.g.
the MK model (for all multi-state partitions) you would type:

```
./raxmlHPC -q part.txt -s picsOrd.phy -K MK -m GTRGAMMA -n T4
```

or to use an ordered model:

```
./raxmlHPC -q part.txt -s picsOrd.phy -K ORDERED -m GTRGAMMA -n T5
```

Note that, it is not possible to assign different multi-state
substitution models to different ambiguous regions, that is, all
multi-state partitions will use the model specified with "-K" or GTR by
default if "-K" is not specified. Moreover, the model of rate
heterogeneity (either CAT or GAMMA) can also not be selected on a per
partition basis, that is, CAT or GAMMA will be used across all
partitions. All remaining RAxML program options (bootstrapping, rapid
bootstrapping, multiple ML searches, estimate of per-partition branch
lengths, etc.) work as for standard RAxML analyses. Thus the only
difference is in the specification of the partition file that is passed
via -q and the new command line switch "-K" that allows you to specify
the desired multi-state substitution model.
