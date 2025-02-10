# Multi-parental QTL analysis

This repository contains the R package mpQTL, which you can find published in [our paper.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04607-z "Multiallelic models for QTL mapping in diverse polyploid populations")

### What can I use `mpQTL` for?

This package can use the [unified mixed model](https://doi.org/10.1038/ng1702 "A unified mixed-model method for association mapping that accounts for multiple levels of relatedness") framework (also known Q + K model) to detect QTLs in genetically diverse populations. Our innovations include support of all polyploid levels and multiallelic markers (like SSRs or haplotypes).

### When is it most useful?

Our research showed that the multiallelic models of `mpQTL` are particularly useful if you have a high number of alleles in your population. For example in multiparental populations, polyploid heterozygous populations or breeding programmes. In that case obtaining multiallelic markers (haplotypes) turned out to be more accurate and powerful than classical biallelic analyses.

## Installation

You can install this package using `devtools::install_github("alethere/mpqtl")`

If you have access to our gitlab and want to install it from there, please contact us for instructions.

## Vignette

You can access the [online vignette here](https://alethere.github.io/mpQTL/ "QTL analysis with the R package mpQTL").

If you want to access your vignette locally you need to build it on install (this will take some time).

```
devtools::install_github("https://github.com/alethere/mpqtl", build_vignette = T)
browseVignettes("mpQTL")
```


