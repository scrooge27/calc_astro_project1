tab=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\tab_gal.txt")
d=tab[,1]
v=tab[,3]

plot(d,v,xlab="distance (Mpc)",ylab="velocity (km/s)")

mod1=lm(v~d)

print(coefficients(mod1))

f1=fitted(mod1)
#fit non-pesato

points(d,f1,col="red",type="l",lwd=2,lty=2)