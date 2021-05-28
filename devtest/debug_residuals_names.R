devtools::load_all()


# import data --------------
str(data)
names(data)
sapply(data, dim)
data$result$phenotype1$beta$PotVar0119975
data$result$phenotype1$residual$PotVar0119912
length(data$result$phenotype1$residual$PotVar0119912)
head(data$pheno)
nrow(data$pheno)


data$snp[1:5,1:5]
class(data$snp)


Kr <- calc.K(t(data$snp))
Kr[1:5,1:5]
head(data$map)




# no NAs in pheno, linear model --------------
test11 <- map.QTL(phenotypes = data$pheno,
                  genotypes = data$snp,
                  ploidy = 4,
                  linear = T,
                  # impute = F, #not needed, since we provided imputed genotypes
                  map = data$map,
                  no_cores = 6)
# There are no residuals in the output (see lm_compare function)
test11$phenotype1$beta$PotVar0119912
test11$phenotype2$beta$PotVar0119912
test11$phenotype2$pval[1:5]


# no NAs in pheno, mixed model --------------
test1 <- map.QTL(phenotypes = data$pheno,
                 genotypes = data$snp,
                 ploidy = 4,
                 K = Kr,
                 K_identity = T, #to force naive
                 # impute = F, #not needed, since we provided imputed genotypes
                 map = data$map,
                 no_cores = 6)
str(head(test1$phenotype1$residual))
str(tail(test1$phenotype1$residual))
test1$phenotype1$residual$PotVar0119912
test1$phenotype2$residual$PotVar0119912
length(test1$phenotype1$residual$PotVar0119912)





# one pheno with NAs, linear model -----------------
phenona <- data$pheno
phenona[1:10,1] <- NA

test22 <- map.QTL(phenotypes = phenona,
                  genotypes = data$snp,
                  ploidy = 4,
                  linear = T,
                  # impute = F, #not needed, since we provided imputed genotypes
                  map = data$map,
                  no_cores = 6)

str(head(test22$phenotype1$pval))
str(head(test22$phenotype1$beta))
str(tail(test22$phenotype1$beta))
str(head(test22$phenotype1$Ftest))
str(head(test22$phenotype1$se))


# one pheno with NAs, mixed model -----------------
test2 <- map.QTL(phenotypes = phenona,
                 genotypes = data$snp[1:9,],
                 ploidy = 4,
                 K = Kr,
                 K_identity = T, #to force naive
                 # impute = F, #not needed, since we provided imputed genotypes
                 map = data$map[1:9,],
                 no_cores = 6)
str(head(test2$phenotype1$residual))
str(tail(test2$phenotype1$residual))
test2$phenotype1$residual$PotVar0119912
test2$phenotype2$residual$PotVar0119912
length(test2$phenotype1$residual$PotVar0119912)
length(test2$phenotype2$residual$PotVar0119912)






####
# phenotypes = data$pheno
# genotypes = data$snp
# ploidy = 4
# K = Kr
# K_identity = T #to force naive
# # impute = F, #not needed, since we provided imputed genotypes
# map = data$map


solveout <- list(
  beta = matrix(-0.5),
  Fstat = 1,
  residual = rep(1,5),
  pval = 1,
  se = 1,
  wald = 1,
  real.df = 1)

solveout[[1]]$residual
result <- list(solveout,solveout,solveout) #example with 3 SNPs
str(result)

# a <- rep(1:nrow(genotypes),each=nrow(phenotypes))
# a[1:500]

#NEW OUTPUT
#Create a list of lists from all results.
#It's actually just 6 lists with k elements where k is markers
res<-do.call(mapply,c(list,result))
res<-as.list(as.data.frame(res))
#All columns are turned into vectors except the beta column
markers <- paste0("snp",1:length(result))

for(r in 2:length(res)) res[[r]] <- unlist(res[[r]])
for(r in 1:length(res)) names(res[[r]]) <- markers

for(r in c(2,4:length(res))) res[[r]] <- unlist(res[[r]])
for(r in c(1:2,4:length(res))) names(res[[r]]) <- markers

res$residual <- split(res$residual,rep(1:3,
                                       each=5))
names(res$residual) <- markers
## for permuted phenotypes store the minimum pvalue only
if (w > npheno) {
  res <- min(res$pval, na.rm = T)
}


###
ordInp_ratios$geno[1:5,1:5]
identical(rownames(ordInp_ratios$geno), ordInp_ratios$map$marker)
identical(colnames(ordInp_ratios$geno), rownames(ordInp_ratios$pheno))
identical(colnames(ordInp_ratios$K), rownames(ordInp_ratios$pheno))


ratios_naive2 <- map.QTL(phenotypes = ordInp_ratios$pheno,
                        genotypes = ordInp_ratios$geno,
                        ploidy = 4,
                        K = ordInp_ratios$K,
                        K_identity = T, #to force naive
                        # impute = F, #not needed, since we provided imputed genotypes
                        map = ordInp_ratios$map,
                        no_cores = 6)

ratios_naive2$pheno01$pval[1:5]
ratios_naive2$pheno01$Fstat[1:9]
ratios_naive2$pheno01$se[1:9]
length(ratios_naive2$pheno01$pval)

names(ratios_naive2$pheno01$residual)[1:9]
length(ratios_naive2$pheno01$residual)
ratios_naive2$pheno01$residual[1:5]
ratios_naive2$pheno01$residual$mrk085432
length(ratios_naive2$pheno01$residual$mrk085432)
