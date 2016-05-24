"""
  @package 
  @file maskSpikes.py
  @mask all the spikes using region files
 
  @brief Created by:
  @author Jun Cheng (Purdue)
 
  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python maskFits .py image.fits region.reg

Notes: 1. 


"""

import pyregion
import numpy as np
import sys
import pyfits
from astropy.io import fits
from astropy.nddata import NDData

def maskFits(regFile, imageFile):
    f1 = open(regFile,'r')
    reg = f1.read() 

    hdulist = pyfits.open(imageFile) 
    h = pyfits.getheader(imageFile)
    hdr = h.copy()
    matrix =hdulist[0].data
    r = pyregion.parse(reg)
    Xrange = matrix.shape[0]
    Yrange = matrix.shape[1]
    
    mask = r.get_mask(hdu=hdulist[0])
    new_matrix = np.zeros((Xrange,Yrange))
    mask_matrix = np.ones((Xrange,Yrange))
    for i in range(Xrange):
        for j in range(Yrange):
            if(mask[i,j]!=True):
                mask_matrix[i,j] = 0 
                new_matrix[i,j]=matrix[i,j]
    # write out the new fits file.
    imageName = imageFile.split('.')
    outName = imageName[0] + "_masked." +imageName[1]
    hdu=fits.PrimaryHDU(new_matrix)
    hdu.writeto(outName,clobber='true')

    hdu_mask=fits.PrimaryHDU(mask_matrix)
    hdu_mask.writeto("mask.fits",clobber='true')

    f1.close()
    hdulist.close()
   
fileName = sys.argv[1]
regFile = sys.argv[2]

maskFits(regFile,fileName)





