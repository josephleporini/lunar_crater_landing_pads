"""Footprint (2.5-D) run: Anderson et al. particle grid from every point of the LM erosion footprint,
each along its own terrain line to the Surveyor 3 camera (vacuum ballistics, first-intersection,
same physics as crater_pads_LIB). Terrain: NAC_DTM_APOLLO12 v1.9 (georef verified vs profile_N49W.csv
to 0.005 m). Slope definitions: 'fwd10' = linear fit of first 10 m toward Surveyor (WP1 definition);
'ctr' = centered plane fit over 5 m radius, directional toward Surveyor."""
import numpy as np, json, sys
sys.path.insert(0,'.')
from geo import *
G=1.625; CAM=1.6; WIN=0.25; Y0=0.9
V=np.arange(10,2390,10.0)
zS=zen(0,0)[0]
eL,nL=en_from_ll(*LM)

def line(e0,n0):
    D=np.hypot(e0,n0); u=np.arange(0,D+0.25,0.25); f=u/D
    z=zen(e0*(1-f),n0*(1-f))-zS          # from source (u=0) to Surveyor (u=D)
    return D,u,z

def slope_fwd10(u,z):
    m=u<=10; return np.rad2deg(np.arctan(np.polyfit(u[m],z[m],1)[0]))

def slope_ctr(e0,n0,r=5.0):
    g=np.arange(-r,r+0.01,0.5); E,N=np.meshgrid(g,g); m=E**2+N**2<=r*r
    Z=zen(e0+E[m],n0+N[m]); a,b,_=np.linalg.lstsq(np.c_[E[m],N[m],np.ones(m.sum())],Z,rcond=None)[0]
    D=np.hypot(e0,n0); return np.rad2deg(np.arctan(a*(-e0/D)+b*(-n0/D)))

def fire(D,u,z,angles):
    T,VV=np.meshgrid(np.deg2rad(angles),V,indexing='ij'); vx=VV*np.cos(T); vy=VV*np.sin(T)
    y0=z[0]+Y0
    tt=u[None,None,:]/vx[...,None]; yp=y0+vy[...,None]*tt-0.5*G*tt**2
    clear=(yp[...,1:-1]>z[None,None,1:-1]).all(-1)     # no ground contact before Surveyor
    yS=yp[...,-1]; hit=clear&(np.abs(yS-(z[-1]+CAM))<=WIN)
    tf=D/vx; vyf=vy-G*tf; sp=np.hypot(vx,vyf); el=np.rad2deg(np.arctan2(vyf,vx))
    return hit, sp, el, T.size

def run(Rs, band, ref):
    pts=[(dn,de) for dn in np.arange(-Rs,Rs+0.1,1.0) for de in np.arange(-Rs,Rs+0.1,1.0) if np.hypot(dn,de)<=Rs]
    tot=0; nh=0; nf=0; pf=0; SP=[]; EL=[]; SL=[]
    for dn,de in pts:
        e0,n0=eL+de,nL+dn; D,u,z=line(e0,n0)
        s={'horizontal':0.0,'fwd10':slope_fwd10(u,z),'ctr':slope_ctr(e0,n0)}[ref]; SL.append(s)
        hit,sp,el,n=fire(D,u,z,band+s); tot+=n; nh+=hit.sum()
        fast=hit&(sp>=100); nf+=fast.sum(); pf+=fast.any(); SP+=list(sp[fast]); EL+=list(el[fast])
    return dict(Rs=Rs, ref=ref, band=[float(band[0]),float(band[-1])], source_points=len(pts),
                slope_median=round(float(np.median(SL)),2) if ref!='horizontal' else 0.0,
                slope_range=[round(float(min(SL)),2),round(float(max(SL)),2)] if ref!='horizontal' else [0,0],
                camera_hits=int(nh), fast_hits=int(nf), points_delivering_fast=int(pf),
                fast_per_1e4=round(1e4*nf/tot,2),
                fast_speed_median=round(float(np.median(SP))) if SP else None,
                fast_arrival_elev_median=round(float(np.median(EL)),2) if EL else None)

if __name__=='__main__':
    D,u,z=line(eL,nL)
    print(f'LM center: fwd10 slope {slope_fwd10(u,z):.2f}, ctr(5 m) slope {slope_ctr(eL,nL):.2f}')
    wide=np.round(np.arange(0,20.01,0.1),2); narrow=np.round(np.arange(1,3.01,0.1),2)
    out={'sweep_center_line_1to3':[], 'footprint':[]}
    # 1) reachability window on the centre line: effective surface slope sweep, 1-3 deg band
    for s in np.arange(-9,-0.99,0.5):
        hit,sp,el,n=fire(D,u,z,narrow+s); f=hit&(sp>=100)
        out['sweep_center_line_1to3'].append([float(s),int(hit.sum()),int(f.sum())])
    print('slope sweep (centre line, 1-3 deg surface band): slope, camera hits, fast hits')
    for r in out['sweep_center_line_1to3']: print('  ',r)
    # 2) footprint runs
    for Rs in (3.0,6.0,10.0):
        for band,ref in [(wide,'horizontal'),(narrow,'horizontal'),(narrow,'fwd10'),(narrow,'ctr'),(wide,'fwd10'),(wide,'ctr')]:
            r=run(Rs,band,ref); out['footprint'].append(r); print(json.dumps(r),flush=True)
    json.dump(out,open('footprint_results.json','w'),indent=1)
