# RERconverge

RERconverge is a set of software written in R that estimates the correlation between relative evolutionary rates of genes and the evolution of a convergent binary or continuous trait across a phylogeny.

## Getting Started

Please refer to the [wiki](https://github.com/nclark-lab/RERconverge/wiki/Getting-Started) for detailed instructions to **install** RERconverge from scratch. For more information on running RERconverge, please see the full documentation ***(Link to R documentation)*** and R vignettes. ***(Link to R vignettes)***

### Quick Start
```
library(devtools)
install_github("nclark-lab/RERconverge")
```
To run an analysis you will need:
1) a trees file: a tab-delimited files with gene names and Newick format trees for each gene.
Tree topologies must be the same for all genes, and at least one tree must contain all species in the dataset.
We provide trees files for several clades [here](bit.ly/2J2QBnj).

2) information about phenotypes for species included in the dataset.
For a *binary trait analysis*, this can either be in the form of:
   * a vector (in R) of species to include in the foreground.
   * a tree object (made in R from a Newick tree) where branches are non-zero only for foreground lineages.

For a *continuous trait analysis*, this should be:
   * a named vector (in R) of quantitative phenotype values, where the names represent the species to which the phenotypes correspond.


### Output

Running RERconverge will produce the following outputs:
1) an object containing, for each gene, the correlation between its relative evolutionary rate and the trait of interest, along with the estimateed p-value
2) an object containing, for each gene, its relative evolutionary rate for each branch of the phylogeny, which can be used in the included visualization scripts (e.g., to illustrate the difference in relative evolutionary rate between foreground and background branches)


## Authors

* **Maria Chikina** - [mchikina](https://github.com/mchikina)
* **Nathan Clark** - [nclark-lab](https://github.com/nclark-lab)
* **Amanda Kowalczyk** - [Amylith](https://github.com/Amylith)
* **Weiguang (Wayne) Mao** - [wgmao](https://github.com/wgmao)
* **Wynn Meyer** - [sorrywm](https://github.com/sorrywm)
* **Raghavendran Partha** - [raghavendranpartha](https://github.com/raghavendranpartha)

See also the list of [contributors](https://github.com/nclark-lab/RERconverge/contributors) who participated in this project.

## Citation

RERconverge can be cited as follows:

```
For coding sequences:
Chikina M, Robinson JD, Clark NL. Hundreds of Genes Experienced Convergent Shifts in Selective Pressure in Marine 
Mammals. Mol Biol Evol. 2016;33: 2182–92. doi:10.1093/molbev/msw112

For conserved non-coding sequences:
Partha R, Chauhan B, Ferreira Z, Robinson J, Lathrop K, Nischal K, et al. Subterranean mammals show convergent 
regression in ocular genes and enhancers, along with adaptation to tunneling. In press eLife. 

```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Methods for computing weights rely on ideas from the following paper:
```
Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read
counts. Genome Biol. 2014;15: R29. doi:10.1186/gb-2014-15-2-r29

```
* Projection operations are drawn from the following paper:
```
Sato T, Yamanishi Y, Kanehisa M, Toh H. The inference of protein-protein interactions by co-evolutionary 
analysis is improved by excluding the information about the phylogenetic relationships. Bioinformatics. 
Bioinformatics Center, Institute for Chemical Research, Kyoto University, Gokasho, Uji, Kyoto 611-0011, 
Japan. sato@kuicr.kyoto-u.ac.jp; 2005;21: 3482–3489. doi:10.1093/bioinformatics/bti564

```
* Thanks to [PurpleBooth](https://github.com/PurpleBooth) for this template.
