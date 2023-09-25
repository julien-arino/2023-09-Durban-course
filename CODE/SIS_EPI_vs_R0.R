# Prevalence vs R_0

# To be able to use LaTeX commands (in the figure labels)
library(latex2exp)

# Values of the EPs
value_EPs = function(R_0, N) {
  EP_I = ifelse(R_0 < 1, 0, (1-1/R_0)*N)
  return(EP_I)
}

R_0 = seq(0.5, 5, by = 0.01)
EP_I = value_EPs(R_0, N = 1000)
# We also show the DFE when R_0>1, so prepare this
R_0_geq_1 = R_0[which(R_0>=1)]
DFE = rep(0, length(R_0_geq_1))

png(file = "../FIGS/endemic_SIS_EE_vs_R0.png",
    width = 1200, height = 800, res = 200)
plot(R_0, EP_I,
     type = "l", lwd = 2,
     xlab = TeX("$R_0$"),
     las = 1,
     ylab = "Prévalence à l'équilibre")
lines(R_0_geq_1, DFE,
      type = "l", lwd = 2,
      lty = 2)
legend("topleft", legend = c("PÉ LAS", "PÉ instable"),
       lty = c(1, 2), lwd = c(2,2),
       bty = "n")
dev.off()
crop_figure(file = "../FIGS/endemic_SIS_EE_vs_R0.png")
