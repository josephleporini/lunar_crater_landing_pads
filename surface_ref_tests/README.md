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
