# coex: An R package for tidy co-expression analyses

---

## How to install

GitHub
```
library(devtools)
devtools::install_github()
```

Bioconductor
```

```

### To do
- Add stopifnot and unit tests for plotPCA
    - Check if plotPCA should be a generic like in DEseq2

- Add spqn normalisation methods. See https://bioconductor.org/packages/release/bioc/vignettes/spqn/inst/doc/spqn.html

- Add CPM/TPM normalisation methods

### Notes

[Measures such as Pearson correlation automatically include a standardization of gene expression across samples, which could explain why additional library size correction may not be needed. This implication is also supported by the Counts workflow outperforming within-sample normalization workflows](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9#Sec9) 

__Signed networks are the default in tidycoex__  
For a nice explanation of why signed networks are usually preferable, see the article [Signed or unsigned: which network type is preferable?](https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/) by Peter Langfelder.

### Further reading
- [Robust normalization and transformation techniques for constructing gene coexpression networks from RNA-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9), [Data + more plots](https://krishnanlab.github.io/RNAseq_coexpression/gtex_plots.html)

- [Addressing confounding artifacts in reconstruction of gene co-expression networks](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1700-9)

- [wTO: weighted topological overlap and a consensus network R package](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2351-7), [Tutorial](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-019-1700-9/MediaObjects/13059_2019_1700_MOESM4_ESM.html)

- Spatial quantile normalisation for co-expression analysis [bioc](https://bioconductor.org/packages/release/bioc/vignettes/spqn/inst/doc/spqn.html#References), [biorxiv](https://www.biorxiv.org/content/10.1101/2020.02.13.944777v1.full.pdf)
