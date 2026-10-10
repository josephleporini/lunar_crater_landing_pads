"""k-sweep cross-check of Topic A Draft 12 Sec. 4.3 using the Anderson et al. (2026) speed grid.
Launch elevation theta = alpha + k*beta (Draft 12 Eq. 2). beta = arctan of least-squares fit over the
first 10 m toward the receptor (Draft 12 definition). Source: 6 m disk, 1 m grid (113 points),
launch height 1.0 m (Draft 12 nominal). Band alpha 1-8 deg and 1-3 deg in 0.1 deg steps.
Speeds 10-2,370 m/s in 10 m/s steps (uniform; Anderson grid), so counts show reachability only.
Hit: within +/-0.25 m of 1.6 m camera height, clear of terrain from 2 m after source to 1 m before receptor.
Usage: A12_DTM=/path/NAC_DTM_APOLLO12.TIF python3 ksweep_run.py"""
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geo import zen, en_from_ll, LM
G=1.62; CAM=1.6; WIN=0.25; Z0=1.0
V=np.arange(10,2380,10.0)
zS=zen(0,0)[0]; eL,nL=en_from_ll(*LM)

def line(e0,n0):
    D=np.hypot(e0,n0); u=np.arange(0,D+0.25,0.25); f=u/D
    return D,u,zen(e0*(1-f),n0*(1-f))-zS

def beta10(u,z):
    m=u<=10; return np.rad2deg(np.arctan(np.polyfit(u[m],z[m],1)[0]))

def fire(D,u,z,theta_deg):
    T,VV=np.meshgrid(np.deg2rad(theta_deg),V,indexing='ij'); vx=VV*np.cos(T); vy=VV*np.sin(T)
    tt=u[None,None,:]/vx[...,None]; yp=z[0]+Z0+vy[...,None]*tt-0.5*G*tt**2
    m=(u>2)&(u<D-1)
    clear=(yp[...,m]>z[None,None,m]).all(-1)
    hit=clear&(np.abs(yp[...,-1]-(z[-1]+CAM))<=WIN)
    vyf=vy-G*(D/vx); sp=np.hypot(vx,vyf)
    return hit, sp

pts=[(dn,de) for dn in np.arange(-6,6.1,1.0) for de in np.arange(-6,6.1,1.0) if np.hypot(dn,de)<=6]
lines=[]; betas=[]
for dn,de in pts:
    D,u,z=line(eL+de,nL+dn); lines.append((D,u,z)); betas.append(beta10(u,z))
out={'source_points':len(pts),'beta_median':round(float(np.median(betas)),2),
     'beta_range':[round(float(min(betas)),2),round(float(max(betas)),2)],'rows':[]}
D0,u0,z0=line(eL,nL); out['beta_touchdown']=round(float(beta10(u0,z0)),2)
for band_name,band in [('1-8',np.round(np.arange(1,8.01,0.1),2)),('1-3',np.round(np.arange(1,3.01,0.1),2))]:
    for k in np.round(np.arange(0,1.21,0.1),1):
        pf=0; nh=0; nf=0; tot=0
        for (D,u,z),b in zip(lines,betas):
            hit,sp=fire(D,u,z,band+k*b); f=hit&(sp>=100)
            pf+=bool(f.any()); nh+=int(hit.sum()); nf+=int(f.sum()); tot+=hit.size
        r=dict(band=band_name,k=float(k),points_fast=pf,hits=nh,fast_hits=nf,fast_share_of_hits=round(nf/nh,3) if nh else None)
        out['rows'].append(r); print(json.dumps(r),flush=True)
print('beta touchdown',out['beta_touchdown'],'disk median',out['beta_median'],'range',out['beta_range'])
json.dump(out,open('ksweep_results.json','w'),indent=1)
