import numpy as np
from scipy.ndimage import map_coordinates
import os, tifffile
A=tifffile.imread(os.environ.get('A12_DTM','NAC_DTM_APOLLO12.TIF')).astype(float); A[A<-1e30]=np.nan
RES=15161.67521207; L0=-44912.5; S0=-2368825.5; CLAT=np.deg2rad(-3.0); R=1737400.0; MPP=2.0000000000006
def rc(lat,lon):   # 0-based row, col (pixel centers)
    return L0 - lat*RES + 1.0, (lon-180.0)*RES*np.cos(CLAT) + S0 + 1.0
def z(lat,lon):
    r,c=rc(np.asarray(lat,float),np.asarray(lon,float))
    return map_coordinates(A,[np.atleast_1d(r),np.atleast_1d(c)],order=1,mode='nearest')
SURV=(-3.0162,336.5820); LM=(-3.0128,336.5781)
# local metric frame centered on Surveyor: east, north (m)
def ll_from_en(e,n):
    return SURV[0]+np.rad2deg(n/R), SURV[1]+np.rad2deg(e/(R*np.cos(np.deg2rad(SURV[0]))))
def en_from_ll(lat,lon):
    return np.deg2rad(lon-SURV[1])*R*np.cos(np.deg2rad(SURV[0])), np.deg2rad(lat-SURV[0])*R
def zen(e,n): return z(*ll_from_en(np.asarray(e,float),np.asarray(n,float)))
