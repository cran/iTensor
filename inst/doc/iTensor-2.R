## ----data, echo=TRUE----------------------------------------------------------
library("iTensor")
data1 <- iTensor::toyModel("MICA")
data2 <- iTensor::toyModel("GroupICA")
str(data1, 2)
str(data2, 2)

## ----data2, echo=TRUE, fig.height=4, fig.width=4------------------------------
plot.ts(data1$X[7700:8000,], main="data1 (X)")
plot.ts(data1$Y[7700:8000,], main="data1 (Y)")

## ----data3, echo=TRUE, fig.height=4, fig.width=8------------------------------
plot.ts(data2$X[4700:5000,], main="data2 (X)")
plot.ts(data2$Y[4700:5000,], main="data2 (Y)")

## ----mica, echo=TRUE, fig.height=4, fig.width=4-------------------------------
t_series <- seq(from = 0.00, to = 1.000, by = 1e-4)
out.MICA <- MICA(data1$X, data1$Y, J=3, gamma_ts = 1 - 1 / (1 + exp(-100 * (t_series - 0.3))))

## ----pairs_mica, echo=TRUE, fig.height=4, fig.width=4-------------------------
plot.ts(out.MICA$U[7700:8000, ], main="Source Signal (X)")
plot.ts(out.MICA$V[7700:8000, ], main="Source Signal (Y)")

## ----groupica, echo=TRUE------------------------------------------------------
out_groupica_pooled_infomax <- GroupICA(data2, J1=6,
    algorithm="pooled", ica.algorithm="InfoMax")

out_groupica_Calhoun2009_fastica <- GroupICA(data2, J1=6,
    algorithm="Calhoun2009", ica.algorithm="FastICA")

out_groupica_Pfister2018_amuse <- GroupICA(data2, J1=6,
    algorithm="Pfister2018", ica.algorithm="AMUSE")

## ----plot_groupica, echo=TRUE, fig.height=4, fig.width=8----------------------
plot.ts(out_groupica_pooled_infomax$Ss[[1]], main="Source Signal (X)")
plot.ts(out_groupica_pooled_infomax$Ss[[2]], main="Source Signal (Y)")

plot.ts(out_groupica_Calhoun2009_fastica$Ss[[1]], main="Source Signal (X)")
plot.ts(out_groupica_Calhoun2009_fastica$Ss[[2]], main="Source Signal (Y)")

plot.ts(out_groupica_Pfister2018_amuse$Ss[[1]], main="Source Signal (X)")
plot.ts(out_groupica_Pfister2018_amuse$Ss[[2]], main="Source Signal (Y)")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

