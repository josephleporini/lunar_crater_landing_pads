import sys, numpy as np
sys.path.insert(0,'/home/claude/lunar_crater_landing_pads')
import crater_pads_LIB as L
# synthetic crater: D=1035 m, d/D=0.15, flat floor radius 38 m, rim at 517.5 m, exterior falls 3 deg
D=1035.; R=D/2; d=0.15*D; rf=38.
x=np.arange(-1600,1601,1.0)
r=np.abs(x)
y=np.where(r<=rf,0.,d*((r-rf)/(R-rf))**1.6)
y=np.where(r>R, d-np.tan(np.deg2rad(3))*(r-R), y)
x=np.concatenate([[-100000],x,[100000]]); y=np.concatenate([[y[0]],y,[y[-1]]])
th=np.rad2deg(np.arctan(np.diff(y)/np.diff(x)))
P=L.Crater_Profile(x,y,th,D,d,'syn',float(th.max()),-R,R)
for x0 in [0,50,-50]:
    print('x0',x0,'signed slope',round(L.signed_surface_slope(x0,P),2))
for vy0 in [5.0,0.0,-5.0]:
    out=[]
    r_=L.run_particle_sim_array([100.0,vy0,50.0,1.0,0.01,10,P,out,D*1.5],mp=False)
    print('vy0',vy0,'returned',type(r_).__name__, 'len' , (len(r_[2]) if isinstance(r_,list) and r_ else 0), 'tracks',len(out))
