ceph=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\ceph_NGC0925c.txt")

tab=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\spline_out_ceph_NGC0925c.txt")

time=ceph[,1]
mag=ceph[,2]

tval=tab[,1]
myspline=tab[,2]

plot(tval,myspline, type='l',xlab="Jd", ylab="mag",col="red",lwd=4)

points(time,mag,col="blue")


