# pip install astroquery matplotlib numpy pandas
import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from astroquery.gaia import Gaia

# Gaia DR3 astrometric binary orbits (Thiele-Innes constants A,B,F,G in mas, period in days).
# The archive caps how many rows one query returns (an Oct 2026 archive update set it to 50,000),
# but ~169k orbits match. So we download in pages of 40,000, sorted by source_id, until a page comes back short.
PAGE = 40000
query = """SELECT TOP {n} source_id, nss_solution_type, period, parallax, a_thiele_innes, b_thiele_innes, f_thiele_innes, g_thiele_innes
           FROM gaiadr3.nss_two_body_orbit
           WHERE a_thiele_innes IS NOT NULL AND parallax / parallax_error > 10 AND source_id > {after}
           ORDER BY source_id"""
os.makedirs('data', exist_ok=True)
path = os.path.join('data', 'gaia_dr3_nss_orbits.csv')
if not os.path.exists(path):                                    # download once (~5 pages, ~1 min each), then reuse
    pages, after = [], -1
    while True:
        page = Gaia.launch_job_async(query.format(n=PAGE, after=after)).get_results().to_pandas()
        pages.append(page); print(f'downloaded page {len(pages)}: {len(page)} rows')
        if len(page) < PAGE: break
        after = page['source_id'].max()
    pd.concat(pages).to_csv(path, index=False)
t = pd.read_csv(path)

A, B, F, G = (t[c + '_thiele_innes'] for c in 'abfg')
u, v = (A**2 + B**2 + F**2 + G**2) / 2, A * G - B * F
a0 = np.sqrt(u + np.sqrt(u**2 - v**2))                          # size of the star's orbital wobble [mas]
t['f_M'] = (a0 / t['parallax'])**3 / (t['period'] / 365.25)**2   # astrometric mass function [Msun]
cand = t[t['f_M'] > 3]               # a DARK companion must be heavier than f_M -> >3 Msun = black-hole territory
print(cand.sort_values('f_M', ascending=False).head(20))
print(t[t['source_id'].isin([4373465352415301632, 5870569352746779008])])   # Gaia BH1 and BH2

plt.scatter(t['period'], t['f_M'], s=1, c='gray', alpha=0.3, label=f'{len(t)} Gaia DR3 orbits')
plt.scatter(cand['period'], cand['f_M'], s=12, c='red', label=f'{len(cand)} with $f_M > 3\\,M_\\odot$')
plt.xscale('log'); plt.yscale('log'); plt.xlabel('Orbital period [days]'); plt.ylabel('Astrometric mass function [$M_\\odot$]')
plt.title('Hunting hidden black holes in Gaia DR3'); plt.legend()
plt.savefig('gaia_dark_companions.png', dpi=150, bbox_inches='tight')
