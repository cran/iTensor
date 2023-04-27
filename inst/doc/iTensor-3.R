## ----data, echo=TRUE----------------------------------------------------------
library("iTensor")
library("nnTensor")
data <- nnTensor::toyModel("CP")
str(data, 2)

## ----data2, echo=TRUE, fig.height=6, fig.width=6------------------------------
plotTensor3D(data)

## ----multilinearica, echo=TRUE------------------------------------------------
out <- MultilinearICA(data, Js=c(4,4,4), algorithm="FastICA")

## ----plot_multilinearica, echo=TRUE, fig.height=6, fig.width=6----------------
rec_data <- recTensor(out$S, out$As)
plotTensor3D(rec_data)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

