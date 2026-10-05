import sys, numpy as np, time
sys.path.insert(0,'/home/claude/lunar_crater_landing_pads')
import crater_pads_LIB as L
D=1035.; R=D/2; d=0.15*D; rf=38.
x=np.arange(-1600,1601,1.0); r=np.abs(x)
y=np.where(r<=rf,0.,d*((r-rf)/(R-rf))**1.6)
y=np.where(r>R, d-np.tan(np.deg2rad(3))*(r-R), y)
x=np.concatenate([[-100000],x,[100000]]); y=np.concatenate([[y[0]],y,[y[-1]]])
th=np.rad2deg(np.arctan(np.diff(y)/np.diff(x)))
P =L.Crater_Profile(x,y,th,D,d,'syn',float(th.max()),-R,R)
Pr=L.Crater_Profile(np.asarray(list(reversed(x)))*-1,np.asarray(list(reversed(y))),np.asarray(list(reversed(th))),D,d,'rev',float(th.max()),-R,R)
v0s=[20,50,100,200,500,1000,2380]; degs=np.arange(1,3.01,0.25)
heights=[1,5,8,16]
def side(profile,x0,surface_ref):
    s=L.signed_surface_slope(x0,profile) if surface_ref else 0.
    bins=np.arange(int(R),int(1.5*D))           # outside the rim on this side
    low=np.full(len(bins),np.inf)
    for v0 in v0s:
        for dg in degs:
            a=np.deg2rad(dg+s); out=[]
            res=L.run_particle_sim_array([v0*np.cos(a),v0*np.sin(a),x0,0.9,0.01,10,profile,out,D*1.5],mp=False)
            if not res: continue
            tx,ty=res[2],res[3]
            m=(tx>=bins[0])&(tx<bins[-1]+1)
            if m.any():
                idx=(tx[m]-bins[0]).astype(int)
                np.minimum.at(low,idx,ty[m])
    ground=np.array([L.calc_crater_wall(b,profile)[0] for b in bins])
    gap=low-ground
    return s,gap
print(f"{'ref':9}{'land x0':>8}{'side':>6}{'slope@x0':>9}{'maxGap(m)':>10}"+''.join(f"{'UZ>'+str(h)+'m':>8}" for h in heights))
t0=time.time()
for ref in [False,True]:
    for x0 in [0,100,200]:
        for name,prof,xx in [('+x',P,x0),('-x',Pr,-x0)]:   # -x side: reversed profile launched from -x0 = same physical point
            s,gap=side(prof,xx,ref)
            g=np.where(np.isinf(gap),999,gap)
            row=f"{'surface' if ref else 'horiz':9}{x0:>8}{name:>6}{s:>9.2f}{min(g.max(),999):>10.1f}"
            row+=''.join(f"{int((g>h).sum()):>8}" for h in heights)
            print(row,flush=True)
print('elapsed s',round(time.time()-t0))
