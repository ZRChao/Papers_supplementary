
## Data generation

- Parameter setting
```r
K <- 30             # the number of the leaf nodes
theta <- 0.1        # the disperison parameters for DM or DM part of ZIDM
pi <- 0.1           # the zero-inflation parameters for ZIDM
beta <- 0.4         # difference between two groups
N1 <- N2 <- 50      # sample size
D1 <- round(runif(N1, 10*K, 1000*K))  # the depth 
D2 <- round(runif(N2, 10*K, 1000*K))
```

- Generate a tree 
```r
library(ape)
set.seed(K)
tree <- rtree(K)
dif_taxa <- 38 # the differential internal nodes
```

- Plot the tree with differential settings
```r
plot(tree, use.edge.length = F, show.tip.label = F, type = "cladogram", direction='downwards')
nodelabels(frame = "circle", cex=0.8, bg='grey')
tiplabels(tip = 1:p, frame = "circle", cex=0.8, bg='grey')

nodelabels(node = dif_taxa, frame = "circle", cex=1, bg='red')
dif_otu <- find_descent(tree$edge, dif_taxa)
(dif_otu <- dif_otu[dif_otu<=p])
tiplabels(tip = dif_otu , frame = "circle", cex=0.8, bg='yellow')
```

- simulate the parameters on the tree------
```r
library(CountOnTree)
para <- set_difpara(tree, theta=theta, dif_taxa=dif_taxa, dif_effect=beta, seed=K, list=F)
head(para)
a1 <- para$para_A
a2 <- para$para_B
```

- Simulate count data on the tree
```r
mnt1 <- rmn_tree(tree, N1, D=D1, para = a1)   # multinomial tree
mnt2 <- rmn_tree(tree, N2, D=D2, para = a2) 
mnt <- rbind(mnt1$otu, mnt2$otu)

dtm_1 <- rdm_tree(tree, N1, D=D1, para = a1)  # Dirichlet tree multinomial
dtm_2 <- rdm_tree(tree, N2, D=D2, para = a2) 
dtm <- rbind(dtm_1$otu, dtm_2$otu)

zidtm_1 <- rzidm_tree(tree, N1, D=D1, para = a1, pi = pi)  # Zero-inflated Dirichlet multinomial 
zidtm_2 <- rzidm_tree(tree, N2, D=D2, para = a2, pi = pi) 
zidtm <- rbind(zidtm_1$otu, zidtm_2$otu)

s <- round(runif(N1+N2, 1, 50))
s[sample(N1+N2, 30)] <- 1
mnt <- round(mnt/s)
dtm <- round(dtm/s)
zidtm <- round(zidtm/s)

group <- c(rep('A', N1), rep('B', N2))
```

## Differential analysis of different methods
```r
dat <- mnt
# t.test for log relative data (0+0.5)
rel_dat <- (dat+0.5)/rowSums(dat+0.5)
tpv <- apply(log(rel_dat), 2, function(x) return(t.test(x~group)$p.value))

# wilcox rank test for raw count data
wpv <- apply(dat, 2, function(x) return(wilcox.test(x~group)$p.value))

# ANCOM 
adat <- data.frame(otu=dat, Group=group)
apv <- ANCOM(adat)$detected

# metagenomeSeq by fitZIG
zpv <- ZIG(dat, group)

# DeSeq2
dpv <- DESEQ2(dat, group)

# edgeR
epv <- EDGER(dat, group)

# adaANCOM
ada <- adaANCOM(tree, dat, group)
summary(ada)
adapv <- ada$otu_info$p.value
```





