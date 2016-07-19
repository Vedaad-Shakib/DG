scal = read.csv("scalability_data.data")
pdf(file="scalability_plot.pdf")
par(mfrow=c(1, 1), lty=2)

# elapsed
lm.fit = lm(elapsed~I(procs^(-1)), data=scal)
plot(procs, elapsed, col="blue")
fitline = predict(lm.fit, data.frame(procs=(.9*min(procs)):(1.1*max(procs))))
lines((.9*min(procs)):(1.1*max(procs)), fitline, col="blue")

elapsed_intercept = signif(coef(summary(lm.fit))["(Intercept)", "Estimate"], digits=4)
elapsed_coef = signif(coef(summary(lm.fit))["I(procs^(-1))", "Estimate"], digits=4)
title(main=paste("Elapsed Time = ", elapsed_intercept, 
                 " + ", elapsed_coef, "procs^(-1)"))

# cpu
lm.fit2 = lm(cpu~I(procs^(-1)), data=scal)
points(procs, cpu, col="green")
fitline = predict(lm.fit2, data.frame(procs=(.9*min(procs)):(1.1*max(procs))))
lines((.9*min(procs)):(1.1*max(procs)), fitline, col="green")

cpu_intercept = signif(coef(summary(lm.fit2))["(Intercept)", "Estimate"], digits=4)
cpu_coef = signif(coef(summary(lm.fit2))["I(procs^(-1))", "Estimate"], digits=4)
title(main=paste("\n\nCPU Time = ", cpu_intercept, 
                 " + ", cpu_coef, "procs^(-1)"))

dev.off()
