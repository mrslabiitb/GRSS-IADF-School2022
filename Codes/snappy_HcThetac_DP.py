# -*- coding: utf-8 -*-
"""
Created on Sun Oct 02 11:09:57 2022

@author: Dr. Dipankar Mandal
References
----------
Mandal, D., Bhattacharya, A., Rao, Y. S., 2021. Radar remote sensing for crop biophysical parameter estimation. Springer.
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
import numpy
import numpy as np
from snappy import Product
from snappy import ProductData
from snappy import ProductIO
from snappy import ProductUtils
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
##---------------------------------------------------------------------------------



##---------------------------------------------------------------------------------
bandc11 = product.getBand('C11')
bandc22 = product.getBand('C22')

##thetac ----------------------------------------------------------------------
thetaProduct = Product('Theta', 'Theta', width, height)
thetaBand = thetaProduct.addBand('theta', ProductData.TYPE_FLOAT32)
HcBand = thetaProduct.addBand('Hc', ProductData.TYPE_FLOAT32)
writer = ProductIO.getProductWriter('BEAM-DIMAP')

ProductUtils.copyGeoCoding(product, thetaProduct)
ProductUtils.copyMetadata(product, thetaProduct)
ProductUtils.copyTiePointGrids(product, thetaProduct)

thetaProduct.setProductWriter(writer)
thetaProduct.writeHeader(name+'_DP_descriptors.dim')





##-----------------------------------------------------------------------------
c11 = numpy.zeros(width, dtype=numpy.float32)
c22 = numpy.zeros(width, dtype=numpy.float32)

print("Writing...")

for y in range(height):
    print("processing line ", y, " of ", height)
    c11 = bandc11.readPixels(0, y, width, 1, c11)
    c22 = bandc22.readPixels(0, y, width, 1, c22)
    
    q = c22/c11
    q[q==1]=1

    mc = (1-q)/(1+q)

    p1 = 1/(1+q)
    p2 = q/(1+q)
    Hc = -1*(p1*np.log2(p1)+p2*np.log2(p2))
    del p1,p2

    thetac = np.arctan(((1-q)**2)/(1-q+q**2)) * (180/np.pi)
    
    #ndvi = (r10 - r7) / (r10 + r7)
    
    theta = thetac
    thetaBand.writePixels(0, y, width, 1, theta)
    HcBand.writePixels(0, y, width, 1, Hc)

    

thetaProduct.closeIO()


print("Done.")