library(rgl)

tab=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\planck_total.txt")
piatti=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\planck_flat3.txt")
chiusi=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\planck_clsd3.txt")
aperti=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_1\\planck_open3.txt")

rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    open3d()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    bg3d(color = bg )
  }
  clear3d(type = c("shapes", "bboxdeco"))
  view3d(theta =0, phi=0, zoom = 0.1)
}


rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "om", ylab="ol", zlab="t", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#E2E2E2","black"))
{ 
  
  lim <- function(x){max(abs(x)) * 1.1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  segments3d(xlim, c(0, 0), c(0, 0), color = axis.col)
  segments3d(c(0, 0), ylim, c(0, 0), color = axis.col)
  segments3d(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  points3d(axes, color = axis.col, size = 3)

text3d(axes, text = c("om", "ol", "res"),
          adj = c(0.5, -0.8), size = 2)

axis3d('x', pos=c( NA, 0, 0 ), col = "darkgrey")
axis3d('y', pos=c( 0, NA, 0 ), col = "darkgrey")
axis3d('z', pos=c( 0, 0, NA ), col = "darkgrey")
}

get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

#generale
om=tab[,1]
ol=tab[,2]
t=tab[,3]

rgl_init()
rgl_add_axes(om, ol, t)

rgl.spheres(om, ol, t, r = 0.0025) 
            #color = get_colors() 

lines3d(-om+1,om,1.25,col="red")

#piatti
om=piatti[,1]
ol=piatti[,2]
t=piatti[,3]

points3d(x=om, y=ol, z=t,col="red")

#chiusi
om=chiusi[,1]
ol=chiusi[,2]
t=chiusi[,3]

points3d(x=om, y=ol, z=t,col="blue")

#aperti
om=aperti[,1]
ol=aperti[,2]
t=aperti[,3]

points3d(x=om, y=ol, z=t,col="green")
