tab=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\chisq_ceph_NGC0925c.txt")

p=tab[,1]
chi=tab[,2]

plot(p,chi,xlab="P (JD)",ylab="chi_square",type="l",lwd=2)

