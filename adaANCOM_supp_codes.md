
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

## Function we used
```r
library(magrittr)
library(MGLM)
library(metagenomeSeq)
ZIG = function(data, group){
    p <- ncol(data)
    if(is.null(colnames(data))) colnames(data) <- as.character(1:p)
    idx<- which(apply(data, 1, function(x) sum(x==0))<(p-1))
    temp_data = t(as.matrix(data[idx, ]))
    group = group[idx]
    indz = AnnotatedDataFrame(data.frame(X = group))
    datz = newMRexperiment(temp_data, phenoData = indz) 
    options(warn=-1)
    qua = cumNormStatFast(datz)
    options(warn=-1)
    new.data = cumNorm(datz, p = qua)#
    
    mod = model.matrix(~1+X, data = pData(new.data))
    options(warn=-1)
    res = fitFeatureModel(new.data, mod)
    res.pv = MRcoefs(res, number = ncol(data))
    
    zig_pv <- res.pv$pvalues[match(colnames(data), rownames(res.pv))]
    return(zig_pv)
}

 library(edgeR)
 EDGER <- function(data, group) {
    temp <- DGEList(counts = t(data), group = factor(group)) # glmQLFit()
    res  <- calcNormFactors(temp) %>% estimateCommonDisp() %>% estimateTagwiseDisp %>% exactTest %>%  topTags(., n=nrow(.))
    res.pv <- res$table$PValue
    res.pv <- res.pv[match(colnames(data), rownames(res$table))]
    res.pv
 }
  
 library(DESeq2)
 DESEQ2 <- function(data, group) {
    warning()
    N <- length(group)
    coldata <- data.frame(row.names = 1:N, condition=factor(group))
    suppressMessages(dds <- DESeqDataSetFromMatrix(countData=t(as.matrix(data)), colData=coldata, design=~condition))
    res <- suppressMessages(DESeq(dds)) %>% results() 
    zig_pv <- res$pvalue[match(colnames(data), rownames(res))]
    return(zig_pv)
 }

 ANCOM <- function(data, group) {
    p <- ncol(data)
    adat <- data.frame(otu=data, Group=group)
    aa <- ANCOM1(adat, theta = 0.1)
    aa <- aa$detected
    colnames(data) -> x
    if(is.null(x)) x <- 1:p
    ancom_pv <- rep(1, p)
    if(length(aa)>0) ancom_pv[match(aa, paste0('otu.', x))] <- 0
    ancom_pv
 }  

ancom.detect1 <- function(otu_data, n_otu, alpha, multcorr, ncore){
  
  ## Detect whether the data are dependent or not
  if( ncol(otu_data) == n_otu+1  ){
    Group     <- otu_data[, ncol(otu_data) ]
    ID        <- rep( 1 , nrow(otu_data) )
    repeated <- FALSE
    fformula  <- formula("lr ~ Group")
  } else if( ncol(otu_data) == n_otu+2  ){
    Group     <- otu_data[, ncol(otu_data)-1 ]
    ID        <- otu_data[, ncol(otu_data)   ]
    repeated <- TRUE
    fformula  <- formula("lr ~ Group | ID")
    
  } else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  ## Detect which test to use, 
  ## Dependent data: Friedman test
  ## Independent data: Wilcoxon Rank-Sum or the Kruskal-Wallis
  ## exactRankTests::wilcox.exact is faster than stats::kruskal.test and stats::wilcox.test
  if( repeated==FALSE ){
    if( length(unique(Group))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  } else{
    tfun <- stats::friedman.test
  }
  
  
  
  ## Parallelized way to get the logratio.mat
  ## Doubles the number of computations to make, so only run the parallel
  ## version if there are multiple cores. Method may also add some computational
  ## overhead, so if only 2 cores, the nested for-loop shoud have advantage
  ## over the parallel loop (though I have not tested that).
  ## For some reason this is taking much longer, do not run the parallel loop as of now.
  if( FALSE ){
    registerDoParallel( cores=ncore )
    
    aa <- bb <- NULL
    logratio.mat <- foreach( bb = 1:n_otu, .combine='rbind', .packages="foreach" ) %:% 
      foreach( aa = 1:n_otu , .combine='c',  .packages="foreach" ) %dopar% {
        if( aa==bb ){
          p_out <- NA
        } else{
          data.pair <- otu_data[,c(aa,bb)]
          lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
          lr_dat <- data.frame( lr=lr, Group=Group, ID=ID )
          p_out  <- tfun(formula=fformula, data = lr_dat)$p.value
        }
        p_out
      }
    rownames(logratio.mat) <- colnames(logratio.mat) <- NULL
  } else{
    logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
    for(ii in 1:(n_otu-1)){
      for(jj in (ii+1):n_otu){
        data.pair <- otu_data[,c(ii,jj)]
        lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
        lr_dat <- data.frame( lr=lr, Group=Group, ID=ID )
        
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }
    }  
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
  }
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==TRUE] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<alpha))
    })
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<alpha))
    })
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<alpha))
    })
  }
  return(W)
  }

############################################################
############################################################




ANCOM1  <-  function(OTUdat, sig=0.05, multcorr=3, tau=0.02, theta=0.1, repeated=FALSE){
  
  #OTUdat <-  read.delim(filepath,header=TRUE)
  
  num_col <- ncol( OTUdat )
  
  if( repeated==FALSE ){
    colnames(OTUdat)[ num_col ] <- "Group"    # rename last column as "Group"
    num_OTU      <- ncol(OTUdat) - 1
    
    sub_drop <- data.frame( nm_drop= "N/A" )
    sub_keep <- data.frame( nm_keep= "All subjects" )
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary   <- paste0( "No subjects entirely removed (not a repeated-measures design)" )
    
  } else{
    colnames(OTUdat)[ num_col-1 ] <- "Group"    # rename 2nd last column as "Group"
    colnames(OTUdat)[ num_col   ] <- "ID"    # rename last column as "ID"
    OTUdat$ID    <- factor( OTUdat$ID )
    num_OTU      <- ncol(OTUdat) - 2
    
    ## Drop subjects if missing at a given time point
    crossTab <- table( OTUdat$Group , OTUdat$ID  )==0
    id_drop  <- apply( crossTab, 2, FUN=function(x) any(x)  )
    nm_drop  <- names( which( id_drop ) )
    idx_drop <- OTUdat$ID %in% nm_drop
    OTUdat   <- OTUdat[ idx_drop==FALSE, ]
    
    if( nrow(OTUdat)==0 ){ stop("Too many missing values in data, all subjects dropped") }    
    OTUdat$ID <- droplevels( OTUdat$ID )    
    num_dropped <- sum(id_drop)
    num_retain  <- length(id_drop) - num_dropped
    
    sub_drop <- data.frame( nm_drop=paste(nm_drop, collapse=", " ) )
    sub_keep <- data.frame( nm_keep= paste(levels(OTUdat$ID), collapse=", " ) )
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary   <- paste0( "Analysis used ", num_retain, " subjects (", num_dropped, " were removed due to incomplete data)")
  }
  
  OTUdat$Group <- factor( OTUdat$Group )
  OTUdat       <- data.frame( OTUdat[ which(is.na(OTUdat$Group)==FALSE),],row.names=NULL )
  
  W.detected   <- ancom.detect1(OTUdat, num_OTU, sig, multcorr, ncore=1 )
  W_stat       <- W.detected
  
  
  
  if( num_OTU < 10 ){
    detected <- colnames(OTUdat)[which(W.detected > num_OTU-1 )]    
  } else{
    ## Detected using a stepwise mode detection
    if( max(W.detected)/num_OTU >= theta ){
      c.start <- max(W.detected)/num_OTU
      cutoff  <- c.start-c(0.05,0.10,0.15,0.20,0.25)
      
      prop_cut <- rep(0,length(cutoff))
      for(cut in 1:length(cutoff)){
        prop_cut[cut] <- length(which(W.detected>=num_OTU*cutoff[cut]))/length(W.detected)
      } 
      
      del <- rep(0,length(cutoff)-1)
      for( ii in 1:(length(cutoff)-1) ){
        del[ii] <- abs(prop_cut[ii]-prop_cut[ii+1])
      }
      
      if(       del[1]< tau & del[2]<tau & del[3]<tau ){ nu=cutoff[1]
      }else if( del[1]>=tau & del[2]<tau & del[3]<tau ){ nu=cutoff[2]
      }else if( del[2]>=tau & del[3]<tau & del[4]<tau ){ nu=cutoff[3]                                
      }else{ nu=cutoff[4] }
      
      up_point <- min(W.detected[ which( W.detected >= nu*num_OTU ) ])
      
      W.detected[W.detected>=up_point] <- 99999
      W.detected[W.detected<up_point]  <- 0
      W.detected[W.detected==99999]    <- 1
      
      detected <- colnames(OTUdat)[which(W.detected==1)]
      
    } else{
      W.detected <- 0
      detected   <- "No significant OTUs detected"
    }
    
  }
  
  #results_list <- list( W         = W_stat,
  #                      Arbitrary = detected_arbitrary,
  #                      Stepwise  = detected_stepwise )
  #idx0 <- lapply( results_list , FUN=length)
  #results_list[idx0==0] <- "No significant OTUs detected"
  
  results <- list( W=W_stat, detected=detected, dframe=OTUdat, repeated=repeated,
                   n_summary=n_summary, sub_drop=sub_drop, sub_keep=sub_keep)
  class(results) <- "ancom"
  
  return(results)  
  
}
```




