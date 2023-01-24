ceph=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\ceph_NGC0925a.txt")

tab=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\spline_out_ceph_NGC0925a.txt")

time=ceph[,1]
app_mag=ceph[,2]
app_mag_err=ceph[,3]


tval=tab[,1]
myspline=tab[,2]

plot(tval,myspline, type='l',xlab="time (JD)", ylab="app_mag",col="red",lwd=3,ylim=c(26.58,24.78))

arrows(x0=time,y0=app_mag-app_mag_err,x1=time,y1=app_mag+app_mag_err,code=3,angle=90,length=0.04,lwd=1.5)

points(time,app_mag,pch=16,cex=0.63)

legend("bottomright",inset=.02,legend="cubic spline",col="red", lty=1,lwd=3)


