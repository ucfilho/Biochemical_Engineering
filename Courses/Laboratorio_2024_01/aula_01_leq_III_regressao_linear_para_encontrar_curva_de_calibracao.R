# substitua seus pontos aqui
x = c(0.1,0.2,0.4,0.5)
y = c(0.2,0.6,1.5,2.0)
xlabel = 'Absorbancia 540nm'
ylabel = 'Concentracao (g/L)'

# Fitting Simple Linear Regression to the Training set
lm.r= lm(y~x)
#Summary of the model
summary(lm.r)

# Create the scatter plot
plot(x, y,
     xlab = xlabel,
     ylab = ylabel,
     main = "Curva de calibração")
abline(lm(y ~ x))
