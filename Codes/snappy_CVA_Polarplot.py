# -*- coding: utf-8 -*-
"""
Created on Sun Oct 02 11:09:57 2022

@author: Dr. Dipankar Mandal

"""

  # ---------------------------------------------------------------------------------------
  # Copyright (C) 2022 by Microwave Remote Sensing Lab, IITBombay http://www.mrslab.in
 
  # This program is free software; you can redistribute it and/or modify it
  # under the terms of the GNU General Public License as published by the Free
  # Software Foundation; either version 3 of the License, or (at your option)
  # any later version.
  # This program is distributed in the hope that it will be useful, but WITHOUT
  # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  # FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
  # more details.
 
  # You should have received a copy of the GNU General Public License along
  # with this program; if not, see http://www.gnu.org/licenses/
  # ---------------------------------------------------------------------------------------
  
import sys
import numpy as np
from snappy import Product
from snappy import ProductData
from snappy import ProductIO
from snappy import ProductUtils
import matplotlib.pyplot as plt
##############


if len(sys.argv) != 2:
    print("usage: %s <file>" % sys.argv[0])
    sys.exit(1)

file = sys.argv[1]


print("Reading...")
product = ProductIO.readProduct(file)
width = product.getSceneRasterWidth()
height = product.getSceneRasterHeight()
name = product.getName()
description = product.getDescription()
band_names = product.getBandNames()

print("Product:     %s, %s" % (name, description))
print("Raster size: %d x %d pixels" % (width, height))
print("Start time:  " + str(product.getStartTime()))
print("End time:    " + str(product.getEndTime()))
print("Bands:       %s" % (list(band_names)))





# Multiplying arrays--------------------------------Plotting

##---------------------------------------------------------------------------------
## https://senbox.atlassian.net/wiki/spaces/SNAP/pages/19300362/How+to+use+the+SNAP+API+from+Python
s1 = product.getBand('rho')
s2 = product.getBand('nu')


# p = ProductIO.readProduct('snappy/testdata/MER_FRS_L1B_SUBSET.dim')
# rad13 = p.getBand('radiance_13')
# w = rad13.getRasterWidth()
# h = rad13.getRasterHeight()
# rad13_data = np.zeros(w * h, np.float32)
# rad13.readPixels(0, 0, w, h, rad13_data)
# p.dispose()
# rad13_data.shape = h, w

w = s1.getRasterWidth()
h = s1.getRasterHeight()
s1_data = np.zeros(w * h, np.float32)
s1.readPixels(0, 0, w, h, s1_data)
s1_data.shape = h, w

s2_data = np.zeros(w * h, np.float32)
s2.readPixels(0, 0, w, h, s2_data)
s2_data.shape = h, w



# ##
xx = s2_data.flatten()
yy = s1_data.flatten()


#histogram definition
bins = [50, 50] # number of bins
# histogram the data
hh, locx, locy = np.histogram2d(xx, yy, bins=bins)
# Sort the points by density, so that the densest points are plotted last
z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(xx,yy)])
idx = z.argsort()
x2, y2, z2 = xx[idx], yy[idx], z[idx]
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
plt.scatter(x2, y2, c=z2, cmap='jet', marker='.')
ax.set_rmax(30)
plt.colorbar()
plt.savefig('CVD_polarplot.png', bbox_inches='tight', dpi=350)

print("Done.")
plt.show()
