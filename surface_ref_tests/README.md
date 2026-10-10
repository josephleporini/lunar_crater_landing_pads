# surface-ref patch: smoke tests

Synthetic 2-D crater (D = 1035 m, d/D = 0.15, flat floor radius 38 m, exterior falls 3 deg).
Not a real DTM. Purpose: verify the patch, not produce results.

Run (no display needed; tkinter is only used for file dialogs):
    mkdir -p stub/tkinter && touch stub/tkinter/__init__.py stub/tkinter/filedialog.py
    PYTHONPATH=stub MPLBACKEND=Agg python3 umbrella_smoke.py

- probe.py          signed slope helper and integrator edge cases (vy0 > 0, = 0, < 0)
- umbrella_smoke.py horizontal vs surface reference, 0-20 deg band, landing offsets (edit x0 list)
- band13.py         same with the measured 1-3 deg band (Immer et al. 2008)

## Surveyor 3 (Apollo 12 LM ejecta)
- surveyor3_run.py: 2-D line profile (profile_N49W.csv), Anderson integrator vs independent check.
- footprint_run.py + geo.py: extended-source (footprint) runs on the full DTM, slope-definition comparison.
  Needs NAC_DTM_APOLLO12.TIF (LROC RDR, not included): A12_DTM=/path/NAC_DTM_APOLLO12.TIF python3 footprint_run.py
  Georeferencing in geo.py verified against profile_N49W.csv to 0.005 m.

## Alignment with Topic A Draft 12 (2026-10-10)
- Vocabulary follows Draft 12: theta_launch = alpha + k*beta. k = 0 is the "horizontal reference",
  k > 0 is "slope-coupled launch", k = 1 is "full slope following". Earlier files here say
  "surface reference" for k = 1.
- Draft 12 beta = forward least-squares fit over the first 10 m toward the receptor
  (crater_pads_LIB.forward_fit_slope; -7.7 deg at the Apollo 12 touchdown point).
- ksweep_run.py reproduces the Draft 12 Sec. 4.3 k thresholds with the Anderson et al. speed grid
  (ksweep_results.json): k <= 0.3 none of 113 points; k 0.4: 9; k 0.5: 63; k 0.6-0.8: 101-110.
- Superseded here: the "centered vs forward-fit slope" framing in footprint_run.py. Draft 12 treats
  the difference as coupling k (about half the slope), not as an ambiguity in beta.
