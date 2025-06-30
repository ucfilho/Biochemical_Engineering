 pH=c(3,4.5,6,7);V=c(64.9,72.3, 84.5, 30.7);
 S=10^(-pH);x1=S;x2=1/S;y=1/V;
 out=nls(y~a+b*x1+c*x2,start=list(a=1,b=2,c=3))
 print(summary(out))
 a = coef(out)[1]; b = coef(out)[2]; c =coef(out)[3]
 Vmax = 1/a
 K1 = 1/b
 K2 = c
 cat("Vmax =", Vmax, ", K1 =", K1, ", K2 =", K2, "\n");
 cat('PH =', log10(K1*K2)/2*(-1), "\n")
