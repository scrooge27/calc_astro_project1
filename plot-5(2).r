tab1=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\planck_flatonly.txt")
tab2=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\planck_flat3.txt")

t0h=1.0911324590290734
tollh=t0h/100

om=tab1[,1]
t=tab1[,3]
plot(om,t,pch=16,cex=0.1,ylim=c(1,1.2))

tollh=tollh*3

polygon(x=c(0,1,1,0),y =c(t0h-tollh, t0h-tollh, t0h+tollh, t0h+tollh),
        col="green", border=NA)
tollh=tollh*2/3
polygon(x=c(0,1,1,0),y =c(t0h-tollh, t0h-tollh, t0h+tollh, t0h+tollh),
        col="blue", border=NA)
tollh=tollh*1/2
polygon(x=c(0,1,1,0),y =c(t0h-tollh, t0h-tollh, t0h+tollh, t0h+tollh),
        col="red", border=NA)

abline(h=t0h, col="black", lty=1,lwd=2.5)


om=tab2[,1]
t=tab2[,3]
points(om,t,pch=21,cex=0.7,col="black",bg="red")