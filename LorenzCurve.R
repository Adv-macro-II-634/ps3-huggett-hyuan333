install.packages("ineq")
library(ineq)
y=Wealth
plot(Lc(y))
ineq(y,type="Gini")
