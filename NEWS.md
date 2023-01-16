# mpQTL 0.6.1
## Bug fixes
- Solved a phenotype naming bug in the mixed model output of **map.QTL**.


# mpQTL 0.6.0

## Minor changes
- In the function **map.QTL** the argument **k** has been renamed **knn**, to
  avoid confusion with the **K** argument used for kinship. 
- **sample.cM** has been renamed **sample.markers** and the **cM** parameter
  renamed **binsize**. 
- **plot.LD** renamed **plot_LD**.

## Bug fixes

- The output residuals table was fixed. It had an incorrect order for
  phenotypes containing missing values.
- Fix a small bug in **plot_LD** causing an error when using less than four 
  quantiles.
- Error in **sample.cM** when using a vector of marker names as input.
- Error in **QQ.plot** using a list of p-values as input.
- Fix bugs in **pheno_box**, that was returning errors when NAs were in
  SNP dosages or haplotypes or phenotypes.
- Fix a bug introduced in a late version of **impute.knn**.
  

## Fulfill CRAN requirements

- Changes in roxygen tags.
- Include examples.
- Resize the example dataset (in data/) to fit in 1 Mb.
- Add documentation for datasets.


# mpQTL 0.3.0
Version distributed to partners at the software workshop in December 2020.

## Minor changes

- Better implementation for continuous genotypes. 
  Changes in functions for input checking and **impute.knn** to allow continuous 
  genotypes.
- More formal input checks. 
  Regarding this, the function **orderInput** has been added, that orders all
  the input matrices. To be used before running **map.QTL**.
- **skyplot** function.
  A few new arguments in the function **skyplot** for more flexibility (e.g. 
  adjust distance between chromosomes, use user-defined colors). Also, 
  plotting time has been considerably reduced.


## Bug fixes

- In the **skyplot** function, a bug related to marker order has been fixed.
- Multiple cofactors passed to the argument **cofactor** of map.QTL were not
  read correctly (only the first one was actually read). 



# mpQTL 0.2.0
Version released for the software workshop in December 2019.
