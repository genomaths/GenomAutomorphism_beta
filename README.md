<!-- README.md is generated from README.Rmd. Please edit that file -->

# GenomAutomorphism [<img src="man/figures/logo.png" align="right" />](https://genomaths.github.io/genomautomorphism)

Robersy Sanchez  
Department of Biology. Eberly College of Science.  
Pennsylvania State University, University Park, PA 16802  
<rus547@psu.edu>  
[ORCID:
orcid.org/0000-0002-5246-1453](https://orcid.org/0000-0002-5246-1453)

## Overview

This is a R package to compute the autimorphisms between pairwise
aligned DNA sequences represented as elements from a Genomic Abelian
group as described in the paper [Genomic Abelian Finite
Groups](https://www.biorxiv.org/content/10.1101/2021.06.01.446543v2). In
a general scenario, whole chromosomes or genomic regions from a
population (from any species or close related species) can be
algebraically represented as a direct sum of cyclic groups or more
specifically Abelian *p*-groups. Basically, we propose the
representation of multiple sequence alignments (MSA) of length *N* as a
finite Abelian group created by the direct sum of homocyclic Abelian
group of *prime-power order*:

   *G* = (ℤ<sub>*p*<sub>1</sub><sup>*α*<sub>1</sub></sup></sub>)<sup>*n*<sub>1</sub></sup> ⊕ (ℤ<sub>*p*<sub>1</sub><sup>*α*<sub>2</sub></sup></sub>)<sup>*n*<sub>2</sub></sup> ⊕ … ⊕ (ℤ<sub>*p*<sub>*k*</sub><sup>*α*<sub>*k*</sub></sup></sub>)<sup>*n*<sub>*k*</sub></sup>

Where, the *p*<sub>*i*</sub>’s are prime numbers, *α*<sub>*i*</sub> ∈ ℕ
and ℤ<sub>*p*<sub>*i*</sub><sup>*α*<sub>*i*</sub></sup></sub> is the
group of integer modulo *p*<sub>*i*</sub><sup>*α*<sub>*i*</sub></sup>.

For the purpose of automorphism between two aligned DNA sequences,
*p*<sub>*i*</sub><sup>*α*<sub>*i*</sub></sup> ∈ {5, 2<sup>6</sup>, 5<sup>3</sup>}.

------------------------------------------------------------------------

## Status

This application is under development. Watch this repo or check for
updates.

------------------------------------------------------------------------

## Tutotials

There are several tutorials on how to use the package at
[GenomAutomorphism](https://genomaths.github.io/genomautomorphism)
website
[<img src="man/figures/logo.png" align="middle" width="32" height="32" />](https://genomaths.github.io/genomautomorphism):

* <a href="https://genomaths.github.io/genomautomorphism/articles/GenomAutomorphism.html" target="_blank">Get started-with GenomAutomorphism</a>
* [Analysis of Automorphisms on a DNA Multiple Sequence Alignment](https://genomaths.github.io/genomautomorphism/articles/automorphism_on_msa.html)
* [Analysis of Automorphisms on a MSA of Primate BRCA1 Gene](https://genomaths.github.io/genomautomorphism/articles/automorphism_on_msa_brca1.html)
* [A Short Introduction to Algebraic Taxonomy on Genes Regions](https://genomaths.github.io/genomautomorphism/articles/automorphism_and_decision_tree.html)
* [Automorphism analysis on COVID-19 data](https://genomaths.github.io/genomautomorphism/articles/covid_19.html)


## Dependences

This package depends, so far, from: *Biostrings*, *GenomicRanges*,
*numbers*, and *S4Vectors*.

------------------------------------------------------------------------

## Installation of R dependencies:

        if (!requireNamespace("BiocManager")) install.packages("BiocManager")
        BiocManager::install()
        
        BiocManager::install(c("Biostrings", "GenomicRanges", "S4Vectors",
        "BiocParallel", "GenomeInfoDb", "BiocGenerics"))
        install.packages(c("numbers", "devtools", "doParallel", "data.table",
        "foreach","parallel"), 
        dependencies=TRUE)

------------------------------------------------------------------------

## You can install **GenomAutomorphism** package from GitHub

       devtools::install_git("https://github.com/genomaths/GenomAutomorphism.git")

------------------------------------------------------------------------

# References

1.  Sanchez R, Morgado E, Grau R. Gene algebra from a genetic code
    algebraic structure. J Math Biol. 2005 Oct;51(4):431-57. doi:
    10.1007/s00285-005-0332-8. Epub 2005 Jul 13. PMID: 16012800. (
    [PDF](https://arxiv.org/pdf/q-bio/0412033.pdf)).

2.  Sanchez R, Grau R, Morgado E. A novel Lie algebra of the genetic
    code over the Galois field of four DNA bases. Math Biosci. 2006;202:
    156–174. <doi:10.1016/j.mbs.2006.03.017>

3.  Sanchez R, Grau R. An algebraic hypothesis about the primeval
    genetic code architecture. Math Biosci. 2009/07/18. 2009;221: 60–76.
    <doi:S0025-5564(09)00114-X> \[pii\] 10.1016/j.mbs.2009.07.001

4.  Robersy Sanchez, Jesús Barreto (2021) Genomic Abelian Finite Groups.
    [doi:
    10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543).

5.  M. V José, E.R. Morgado, R. Sánchez, T. Govezensky, The 24 possible
    algebraic representations of the standard genetic code in six or in
    three dimensions, Adv. Stud. Biol. 4 (2012)
    119–152.[PDF](https://is.gd/na9eap).

6.  R. Sanchez. Symmetric Group of the Genetic–Code Cubes. Effect of the
    Genetic–Code Architecture on the Evolutionary Process MATCH Commun.
    Math. Comput. Chem. 79 (2018) 527-560.
    [PDF](https://bit.ly/2Z9mjM7).

## See also

[Symmetric Group of the Genetic-Code
Cubes](https://github.com/genomaths/GenomeAlgebra_SymmetricGroup)
