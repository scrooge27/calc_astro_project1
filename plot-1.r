#leggo il file originale e fitto la retta-STEP 1
ceph=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\ceph_catalog.txt")

p=ceph[,1]
abs_mag=ceph[,2]

p=log10(p)

mod1=lm(abs_mag~p)

print(coefficients(mod1))

f1=fitted(mod1)

plot(p,f1,col="blue",type="l",lwd=2,lty=1,xlab="log(P)",ylab="abs_mag",ylim=c(-1,-6))

points(p,abs_mag,pch=16,cex=1.3,col="red")



#leggo il file generato e plotto i punti-STEP 3
mags=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\tab_ceph.txt")

calc_p=mags[,1]
calc_p_err=mags[,2]
calc_m=mags[,5]
calc_m_err=mags[,6]

calc_p=log10(calc_p)

points(calc_p,calc_m,pch=16,cex=1.1)
#arrows(x0=calc_p,y0=calc_m-calc_m_err,x1=calc_p,y1=calc_m+calc_m_err,code=3,angle=90,length=0.05,lwd=1.5)
arrows(x0=calc_p-calc_p_err,y0=calc_m,x1=calc_p+calc_p_err,y1=calc_m,code=3,angle=90,length=0.05,lwd=1.5)

legend("bottomright",inset=.05,legend=c("linear fit","cefeidi fornite","cefeidi calcolate"),col=c("blue","red","black"), lty=1,pch=c(20,16,16),lwd=1,cex=0.8)
