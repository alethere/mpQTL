
# mpQTL 0.5.0

## Bug fixes

- The displacement of residuals was fixed. It was incorrect when a phenotype
  contained missing values.
- Error in **sample.cM** when using a vector of marker names as input.
- Error in **QQ.plot** using a list of p-values as input.
- Fix bugs in **pheno_box**, that was returning errors when NAs were in
  SNP dosages or haplotypes or phenotypes.

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
