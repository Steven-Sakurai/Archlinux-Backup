library(ISLR)
attach(Wage)

## poly regression

fit = lm(wage ~ poly(age, 4, raw = T), data = Wage)
coef(summary(fit))

## regression splines

library(splines)
fit1 = lm(wage ~ bs(age, knots = c(25, 40, 60)), data = Wage)
fit2 = lm(wage ~ bs(age, df = 4, degree = 3), data = Wage)
agelim = range(age)
age.grid = seq(agelim[1], agelim[2])
pred1 = predict(fit, newdata = list(age = age.grid), se = T)
pred2 = predict(fit2, newdata = list(age = age.grid), se = T)

setEPS()
postscript("test.eps")

plot(age, wage, col = "grey")
lines(age.grid, pred1$fit, lwd = 2)
lines(age.grid, pred1$fit + 2*pred1$se, lty = "dashed")
lines(age.grid, pred1$fit - 2*pred1$se, lty = "dashed")
lines(age.grid, pred2$fit, lwd = 2, col = "blue")
lines(age.grid, pred2$fit + 2*pred2$se, lty = "dashed", col = "blue")
lines(age.grid, pred2$fit - 2*pred2$se, lty = "dashed", col = "blue")

dev.off()





