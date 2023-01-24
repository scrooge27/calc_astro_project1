import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import loadtxt

os.chdir("C://Users//simop//desktop//uni//ii anno//calc-astro//project_1")
x,y,z= loadtxt('planck.txt',unpack=True)


sc=plt.scatter(x, y, marker='+', s=150, linewidths=4, c=y, cmap="viridis")
cb=plt.colorbar(sc)

plt.show()