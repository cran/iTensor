## ----data, echo=TRUE----------------------------------------------------------
library("iTensor")
data1 <- iTensor::toyModel("ICA_Type1")
data2 <- iTensor::toyModel("ICA_Type2")
data3 <- iTensor::toyModel("ICA_Type3")
data4 <- iTensor::toyModel("ICA_Type4")
data5 <- iTensor::toyModel("ICA_Type5")
dim(data1$X_observed)
dim(data2$X_observed)
dim(data3$X_observed)
dim(data4$X_observed)
dim(data5$gene) # N < P systems

## ----data2, echo=TRUE, fig.height=4, fig.width=4------------------------------
pairs(data1$X_observed, main="data1")
pairs(data2$X_observed, main="data2")
pairs(data3$X_observed, main="data3")
pairs(data4$X_observed, main="data4")
dim(data5$gene)

## ----setting, echo=TRUE-------------------------------------------------------
allICA <- function(X, J){
    # Classic ICA Methods
    out.FastICA <- ICA(X=X, J=J, algorithm="FastICA")
    out.InfoMax <- ICA(X=X, J=J, algorithm="InfoMax")
    out.ExtInfoMax <- ICA(X=X, J=J, algorithm="ExtInfoMax")
    # Modern ICA Method
    out.JADE <- ICA2(X=X, J=J, algorithm="JADE")
    # Output
    list(out.FastICA=out.FastICA, out.InfoMax=out.InfoMax,
      out.ExtInfoMax=out.ExtInfoMax, out.JADE=out.JADE)
}

CorrIndexAllICA <- function(S, out){
    tmp <- lapply(out, function(x, S){CorrIndex(cor(S, Re(x$S)))}, S=S)
    Name <- gsub("out.", "", names(tmp))
    Value <- round(unlist(tmp), 3)
    names(Value) <- NULL
    cbind(Name, Value)
    tmp <- as.data.frame(cbind(Name, Value))
    colnames(tmp) <- c("Method", "CorrIndex")
    knitr::kable(tmp)
}

## ----ica1, echo=TRUE----------------------------------------------------------
set.seed(123456)
out1 <- allICA(X=data1$X_observed, J=3)

## ----corrindex1, echo=TRUE----------------------------------------------------
CorrIndexAllICA(data1$S_true, out1)

## ----pairs_ica1, echo=TRUE, fig.height=4, fig.width=4-------------------------
pairs(out1$out.FastICA$S, main="FastICA")
pairs(out1$out.InfoMax$S, main="InfoMax")
pairs(out1$out.ExtInfoMax$S, main="ExtInfoMax")
pairs(Re(out1$out.JADE$S), main="JADE")

## ----ica2, echo=TRUE----------------------------------------------------------
set.seed(123456)
out2 <- allICA(X=data2$X_observed, J=3)

## ----corrindex2, echo=TRUE----------------------------------------------------
CorrIndexAllICA(data2$S_true, out2)

## ----pairs_ica2, echo=TRUE, fig.height=4, fig.width=4-------------------------
pairs(out2$out.FastICA$S, main="FastICA")
pairs(out2$out.InfoMax$S, main="InfoMax")
pairs(out2$out.ExtInfoMax$S, main="ExtInfoMax")
pairs(Re(out2$out.JADE$S), main="JADE")

## ----ica3, echo=TRUE----------------------------------------------------------
set.seed(123456)
out3 <- allICA(X=data3$X_observed, J=3)

## ----corrindex3, echo=TRUE----------------------------------------------------
CorrIndexAllICA(data3$S_true, out3)

## ----pairs_ica3, echo=TRUE, fig.height=4, fig.width=4-------------------------
pairs(out3$out.FastICA$S, main="FastICA")
pairs(out3$out.InfoMax$S, main="InfoMax")
pairs(out3$out.ExtInfoMax$S, main="ExtInfoMax")
pairs(Re(out3$out.JADE$S), main="JADE")

## ----ica4, echo=TRUE----------------------------------------------------------
set.seed(123456)
out4 <- allICA(X=data4$X_observed, J=3)

## ----corrindex4, echo=TRUE----------------------------------------------------
CorrIndexAllICA(data4$S_true, out4)

## ----pairs_ica4, echo=TRUE, fig.height=4, fig.width=4-------------------------
pairs(out4$out.FastICA$S, main="FastICA")
pairs(out4$out.InfoMax$S, main="InfoMax")
pairs(out4$out.ExtInfoMax$S, main="ExtInfoMax")
pairs(Re(out4$out.JADE$S), main="JADE")

## ----ica5, echo=TRUE----------------------------------------------------------
library("mixOmics")
set.seed(123456)

# mixOmics's IPCA
res.ipca <- ipca(data5$gene, ncomp=3, mode="deflation")

# iTensor's IPCA
out5 <- ICA2(X=as.matrix(data5$gene), J=3, algorithm="IPCA")

## ----pairs_ica5, echo=TRUE, fig.height=4, fig.width=4-------------------------
pairs(res.ipca$x, main="IPCA (mixOmics)", col=data5$treatment[,"Treatment.Group"], pch=16)
pairs(out5$S, main="IPCA (iTensor)", col=data5$treatment[,"Treatment.Group"], pch=16)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

