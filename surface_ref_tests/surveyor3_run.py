#!/usr/bin/env python3
"""Apollo 12 LM -> Surveyor 3: Anderson et al. (2026) ballistics on the measured terrain.

Terrain: profile_N49W.csv (LROC NAC_DTM_APOLLO12 v1.9, 2 m/px, bilinear at 1 m), along
N 48.9 W from Surveyor 3 (dist 0) to the Apollo 12 LM (dist 156.8 m). Built from the
LROC-measured positions (Wagner et al. 2017). See Surveyor 3 Validation Data Sheet rev. 7.

Frame used here: x' = 156.8 - dist, so the LM is at x' = 0 and Surveyor at x' = 156.8,
and LM ejecta toward Surveyor travel in +x'. Heights relative to DTM ground at Surveyor.

Two integrators:
  A) crater_pads_LIB.run_particle_sim_array (Anderson et al., with the surface-ref patch's
     downward-launch guard). Its time window can end a track in mid-air; see check B.
  B) Independent check: same vacuum ballistics (g = 1.625), first ground intersection found
     on a 0.25 m x-grid out to Surveyor. No time window.
Receptor: Surveyor 3 TV camera, 1.6 m (+/- 0.25 m) above ground at x' = 156.8.

Usage: PYTHONPATH=<tkinter stub> MPLBACKEND=Agg python3 surveyor3_run.py profile_N49W.csv
"""
import sys, os, numpy as np, multiprocessing as mp, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import crater_pads_LIB as L

G = L.g
SEP = 156.8           # LM to Surveyor [m]
CAM = 1.6             # camera height above Surveyor ground [m]
CAM_TOL = 0.25
Y0 = 0.9              # launch height above local ground (Anderson et al. Methods)

d = np.loadtxt(sys.argv[1] if len(sys.argv) > 1 else 'profile_N49W.csv', delimiter=',', comments='#')
dist, h = d[:, 0], d[:, 1]
xp = SEP - dist
order = np.argsort(xp)
xp, h = xp[order], h[order]
# pad the ends the way crater_pads_LIB does (flat extension to +/-100 km)
X = np.concatenate([[-100000.], xp, [100000.]])
Y = np.concatenate([[h[0]], h, [h[-1]]])
TH = np.rad2deg(np.arctan(np.diff(Y) / np.diff(X)))
WIDTH = 2 * (xp.max())          # puts the code's "wall" check at the far end of the profile
PROF = L.Crater_Profile(X, Y, TH, WIDTH, float(h.max() - h.min()), 'A12_N49W', float(np.abs(TH).max()), xp.min(), xp.max())

def ground(x):
    return np.interp(x, xp, h)

GROUND_SURV = float(ground(SEP))
CAM_LO, CAM_HI = GROUND_SURV + CAM - CAM_TOL, GROUND_SURV + CAM + CAM_TOL

def slope_toward_surveyor(half_window):
    return float(np.rad2deg(np.arctan((ground(half_window) - ground(-half_window)) / (2 * half_window))))

# ---------- integrator A: Anderson et al. code ----------
def run_A(args):
    v0, ang = args
    a = np.deg2rad(ang)
    res = L.run_particle_sim_array([v0*np.cos(a), v0*np.sin(a), 0.0, Y0, 0.01, 10, PROF, [], WIDTH*1.5], mp=False)
    if not res:
        return (v0, ang, 'dropped', np.nan, np.nan, np.nan)
    vx, vy, x, y = res[0], res[1], res[2], res[3]
    if x[-1] < SEP:                                  # track ended before Surveyor
        gx = ground(x[-1])
        status = 'impact' if y[-1] - gx < 1.0 else 'truncated'
        return (v0, ang, status, float(x[-1]), np.nan, np.nan)
    yS = float(np.interp(SEP, x, y))
    i = int(np.searchsorted(x, SEP)); i = min(i, len(vx)-1)
    sp = float(np.hypot(vx[i], vy[i])); el = float(np.rad2deg(np.arctan2(vy[i], vx[i])))
    return (v0, ang, 'crosses', yS, sp, el)

# ---------- integrator B: independent check ----------
XG = np.arange(0.0, SEP + 0.25, 0.25)
HG = ground(XG)
def run_B(args):
    v0, ang = args
    a = np.deg2rad(ang)
    vx, vy = v0*np.cos(a), v0*np.sin(a)
    y0 = ground(0.0) + Y0
    t = XG / vx
    y = y0 + vy*t - 0.5*G*t**2
    hit = np.where(y <= HG)[0]
    if len(hit):
        return (v0, ang, 'impact', float(XG[hit[0]]), np.nan, np.nan)
    yS = float(y[-1]); vyS = vy - G*t[-1]
    return (v0, ang, 'crosses', yS, float(np.hypot(vx, vyS)), float(np.rad2deg(np.arctan2(vyS, vx))))

def summarize(rows, label):
    cross = [r for r in rows if r[2] == 'crosses']
    hits_cam = [r for r in cross if CAM_LO <= r[3] <= CAM_HI]
    hits_body = [r for r in cross if GROUND_SURV <= r[3] <= GROUND_SURV + 3.0]
    gap = (min(r[3] for r in cross) - GROUND_SURV) if cross else float('inf')
    out = {
        'case': label, 'n': len(rows),
        'dropped': sum(r[2] == 'dropped' for r in rows),
        'truncated_midair': sum(r[2] == 'truncated' for r in rows),
        'impact_before_surveyor': sum(r[2] == 'impact' for r in rows),
        'cross_surveyor': len(cross),
        'lowest_ejecta_above_surveyor_ground_m': round(gap, 2),
        'hit_camera_window': len(hits_cam),
        'hit_spacecraft_0_to_3m': len(hits_body),
    }
    if hits_cam:
        sp = np.array([r[4] for r in hits_cam]); el = np.array([r[5] for r in hits_cam])
        ang = np.array([r[1] for r in hits_cam])
        out.update({'cam_speed_min': round(sp.min()), 'cam_speed_median': round(float(np.median(sp))),
                    'cam_frac_speed_ge_100': round(float((sp >= 100).mean()), 3),
                    'cam_arrival_elev_deg_median': round(float(np.median(el)), 2),
                    'cam_launch_angle_range_deg': [round(ang.min(), 2), round(ang.max(), 2)]})
    return out

if __name__ == '__main__':
    t0 = time.time()
    v0s = np.arange(10, 2380 + 10, 10, dtype=float)          # Anderson et al. grid
    wide = np.round(np.arange(0, 20 + 0.1, 0.1), 2)            # Anderson et al. grid
    narrow = np.round(np.arange(1, 3 + 0.1, 0.1), 2)           # Immer et al. 2008
    s_line = slope_toward_surveyor(5.0)                        # 10 m central difference along N49W
    s_wp1 = -7.0                                               # WP1 footprint median (-6.9 to -7.6)
    print(f'Ground at LM {ground(0):.2f} m, at Surveyor {GROUND_SURV:.2f} m; camera window {CAM_LO:.2f}-{CAM_HI:.2f} m')
    print(f'Slope toward Surveyor at LM: line 10 m fit {s_line:.2f} deg; WP1 footprint {s_wp1:.1f} deg')
    cases = [
        ('horizontal 0-20 (Anderson et al.)',       wide,   0.0),
        ('surface 0-20, line slope',                wide,   s_line),
        ('surface 0-20, WP1 slope -7',              wide,   s_wp1),
        ('horizontal 1-3 (Immer band)',             narrow, 0.0),
        ('surface 1-3, line slope',                 narrow, s_line),
        ('surface 1-3, WP1 slope -7',               narrow, s_wp1),
    ]
    results = []
    with mp.Pool(max(1, os.cpu_count())) as pool:
        for label, band, s in cases:
            args = [(v, round(a + s, 4)) for v in v0s for a in band]
            rA = pool.map(run_A, args, chunksize=200)
            rB = pool.map(run_B, args, chunksize=500)
            sa, sb = summarize(rA, label + ' | A: Anderson code'), summarize(rB, label + ' | B: check')
            results += [sa, sb]
            for s_ in (sa, sb):
                print(json.dumps(s_), flush=True)
    with open('surveyor3_results.json', 'w') as f:
        json.dump({'slope_line_deg': s_line, 'slope_wp1_deg': s_wp1, 'results': results}, f, indent=1)
    print('elapsed s', round(time.time() - t0))
