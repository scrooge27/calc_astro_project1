tab=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\tab_gal.txt")
dp=tab[,1]
d_err=tab[,2]
v=tab[,3]
v_err=tab[,4]
h=tab[,5]
h_err=tab[,6]

d=c(0:26)

plot(d,d*h[1],type="l",lwd=0.5,lty=5,col="blue",xlab="distance (Mpc)",ylab="velocity (km/s)",xlim=c(0,25),ylim=c(0,2000))
points(d,d*h[2],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[3],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[4],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[5],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[6],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[7],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[8],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[9],type="l",lwd=0.5,lty=5,col="blue")
points(d,d*h[10],type="l",lwd=0.5,lty=5,col="blue")



#la retta media
h0=77.37736
points(d,d*h0,type="l",col="red",lwd=2,lty=1)

points(dp,v,pch=16)
arrows(x0=dp,y0=v-v_err,x1=dp,y1=v+v_err,code=3,angle=90,length=0.07,lwd=1.5)
arrows(x0=dp-d_err,y0=v,x1=dp+d_err,y1=v,code=3,angle=90,length=0.07,lwd=1.5)

legend("bottomright",inset=.05,legend=c("galaxies","mean"), lty=c(5,1),lwd=c(1,3),col=c("blue","red"))



#la retta teorica
#mod1=lm(v~dp)

#print(coefficients(mod1))

#f1=fitted(mod1)
#fit non-pesato

#points(dp,f1,col="red",type="l",lwd=2,lty=2)

text(17,2000,"101",font=2)
text(25,1450,"63",font=2)
text(23,1900,"77",font=2)