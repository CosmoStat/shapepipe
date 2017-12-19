#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:04:47 2017

@author: aguinot
"""

"""!
    COMMAND LINE

    python create_flags2.py ../CFIS.093.301.r-000-0.fits --halo_maskmodel halo_mask.reg --spike_maskmodel MEGAPRIME_star_i_13.8.reg

"""




import argparse
import sys
import subprocess
import astropy.coordinates as coord
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import numpy as np
import re
import os

#import matplotlib.pyplot as plt


def sphereDist(p1, p2, wcs):
    """!
        Compute spheric distance between 2 points
        param p1 array [x,y] first point
        param p2 array [x,y] second point
        param wcs astropy.wcs object containing wcs image informations
        return the distance in degree
    """
    p1 = (np.pi/180.)*np.hstack(wcs.all_pix2world (p1[0], p1[1], 1))
    p2 = (np.pi/180.)*np.hstack(wcs.all_pix2world (p2[0], p2[1], 1))

    dTheta = p1 - p2
    dLat = dTheta[1] #dec
    dLong = dTheta[0] #ra

    dist = 2*np.arcsin(np.sqrt(np.sin(dLat/2)**2 +
                               np.cos(p1[1])*np.cos(p2[1])*np.sin(dLong/2)**2))

    return dist*(180./np.pi)*3600.


def find_stars(ra, dec, r):
    """!
        Return GSC (Guide Star Catalog) objects for a field with center (ra,dec) and radius r
        param ra right ascention astropy.wcs oject
        param dec declinaison astropy.wcs oject
        param r radius in arcmin
        return stars dicotionnary for GSC objects in the field
    """
    if dec>0:
        sign='+'
    else:
        sign=''

    s=subprocess.check_output('findgsc2.2 {0} {1}{2} -r {3} -n 1000000'.format(ra, sign, dec, r), shell=True)

    h=[]
    stars={}
    #get 'header'
    for i in s.splitlines()[3].split(' '):
        if (i != '') & (i!=';'):
            #clean output
            i=i.replace(' ','')
            for v in re.split(',|#|;',i):
                if v!='':
                    i=v

            h.append(i)
            stars[i]=[]

    #get 'data'
    for i in range(4,len(s.splitlines())-5):
        k=0
        for j in s.splitlines()[i].split(' '):
            if (j != '') & (j!=';'):
                #clean output
                j=j.replace(' ','')
                for v in re.split(',|#|;',j):
                    if v!='':
                        j=v
                #handle missing data
                try:
                    j=float(j)
                    stars[h[k]].append(j)
                except:
                    if j=='---':
                        stars[h[k]].append(None)
                    else:
                        stars[h[k]].append(j)
                k+=1

    return stars


def apply_mask(input_maskmodel_path, wcs, stars, mag_limit, output_mask_path, mag_pivot=13.8, scale_factor=0.3):
    """!
        Apply mask from model to stars and save into DS9 region file
        @param input_maskmodel_path path to mask model for spike
        @param wcs astropy.wcs object containing wcs image informations
        @param stars stars dico (output of find_stars)
        @param mag_limit higher magnitude to apply the mask
        @param output_mask_path path to save the output region file
        @param mag_pivot pivot magnitude for the model
        @param scale_factor scaling for the model
    """
    mask_file=np.loadtxt(input_maskmodel_path).transpose()
    mask=open(output_mask_path,'w')

    stars_used=[[],[],[]]
    for ra,dec,Fmag, Jmag,Vmag,Nmag,clas in zip(stars['RA(J2000)'],stars['Dec(J2000)'],stars['Fmag'], stars['Jmag'], stars['Vmag'], stars['Nmag'], stars['Clas']):
        mag=0.
        i=0.
        if Fmag!=None:
            mag+=Fmag
            i+=1.
        if Jmag!=None:
            mag+=Jmag
            i+=1.
        if Vmag!=None:
            mag+=Vmag
            i+=1.
        if Nmag!=None:
            mag+=Nmag
            i+=1.
        if i==0.:
            mag=None
        else:
            mag/=i

        if (ra!=None) & (dec!=None) & (mag!=None) & (clas!=None):
            if (mag<mag_limit) & (clas==0):
                scaling = 1. - scale_factor * (mag - mag_pivot)
                pos = wcs.all_world2pix(ra,dec,0)
                stars_used[0].append(pos[0])
                stars_used[1].append(pos[1])
                stars_used[2].append(scaling)



    for i in range(len(stars_used[0])):
        poly = 'polygon('
        for x,y in zip(mask_file[0],mask_file[1]):
            angle = np.arctan2(y,x)
            l = stars_used[2][i] * np.sqrt(x * x + y * y)
            xnew = l * np.cos(angle)
            ynew = l * np.sin(angle)

            poly = poly + str(stars_used[0][i] + xnew + 0.5) + ' ' + str(stars_used[1][i] + ynew + 0.5) + ' '
        poly = poly + ')\n'
        mask.write(poly)

    mask.close()



def main(argv):

    parser = argparse.ArgumentParser(description='Create a flag_map for CFIS images')
    parser.add_argument('imagefile', help='Image for which masks should be made')
    parser.add_argument('--brightstar-maglimit', type=float, help='Faintest stars to mask with complicated mask',
                        default = 18.0, dest='brightstar_maglimit')
    parser.add_argument('--halo_maskmodel', help='File containting halo masking model')
    parser.add_argument('--spike_maskmodel', help='File containting spike masking model')
#    parser.add_argument('outputfile', help = 'Name of output flag_map.')

    args = parser.parse_args()


    img = pyfits.open(args.imagefile)
    header = img[0].header
    img.close()
    wcs=pywcs.WCS(header)



    fieldcenter_pix= np.array([header['CRPIX1'], header['CRPIX2']])
    fieldcenter_wcs= coord.SkyCoord(ra=header['CRVAL1'], dec=header['CRVAL2'], unit='deg')

    radius = sphereDist(fieldcenter_pix, np.zeros(2), wcs)/3600.

    stars=find_stars(fieldcenter_wcs.ra.value, fieldcenter_wcs.dec.value, radius*60.)

#    spikemask_path='MEGAPRIME_star_i_13.8.reg'
    spikemask_path=args.spike_maskmodel
#    halomask_path='halo_mask.reg'
    halomask_path=args.halo_maskmodel

    output_spikemask_path='spike.reg'
    output_halomask_path='halo.reg'

    #-p 12. -m 0.1 -l 18.0

    #Mask spikes
    apply_mask(spikemask_path, wcs, stars, 18., output_spikemask_path)

    #Mask halos
    apply_mask(halomask_path, wcs, stars, 13., output_halomask_path, scale_factor=0.05)

    flag=np.zeros((int(wcs.to_header()['CRPIX1']*2),int(wcs.to_header()['CRPIX2']*2)),dtype=int)

    width=300
    flag[0:width,:]=4
    flag[-width:-1,:]=4
    flag[:,0:width]=4
    flag[:,-width:-1]=4

    #f=pyfits.ImageHDU(flag)
    #f.writeto('test_flag.fits',overwrite=True)
#    weightfile='{0}.{1}.fits'.format(re.split("\-(.*[0-9])\-(.*[0-9])\.",args.imagefile)[0],'weight')
    weightfile='{0}.{1}.fits'.format(re.split(".fits",args.imagefile)[0],'weight')

    os.system('ww -c default.ww -WEIGHT_NAMES {0} -POLY_NAMES halo.reg -POLY_OUTFLAGS 1 -FLAG_NAMES "" -OUTFLAG_NAME flag_CFIS.fits -OUTWEIGHT_NAME ""'.format(weightfile))
    os.system('ww -c default.ww -WEIGHT_NAMES {0} -POLY_NAMES spike.reg -POLY_OUTFLAGS 2 -FLAG_NAMES "" -OUTFLAG_NAME flag_CFIS2.fits -OUTWEIGHT_NAME ""'.format(weightfile))

    hf=pyfits.open('flag_CFIS.fits')
    sf=pyfits.open('flag_CFIS2.fits')

    final_flag=flag+hf[0].data[:,:]+sf[0].data[:,:]
    final_flag=final_flag.astype(np.int16,copy=False)
    f=pyfits.PrimaryHDU(final_flag)
    os.system('rm flag_CFIS.fits')
    os.system('rm flag_CFIS2.fits')

#    flagfile='{0}.{1}.fits'.format(re.split("\-(.*[0-9])\-(.*[0-9])\.",args.imagefile)[0],'flag')
    flagfile='{0}.{1}.fits'.format(re.split(".fits",args.imagefile)[0],'flag')
    f.writeto(flagfile,overwrite=True)

##########

#if __name__ == '__main__':
#    main(argv = sys.argv)
