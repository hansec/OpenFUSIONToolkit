"""Compare railgun plot-file fields (same mesh/decomposition) against the baseline at matching times."""
import glob, sys, os
import numpy as np
import h5py


def load(d):
    files = sorted(glob.glob(os.path.join(d, 'oft_xdmf.*.h5')))
    out = {}
    for f in files:
        with h5py.File(f, 'r') as h:
            m = h['mug/vmesh']
            for k in m.keys():
                if not k.isdigit() or int(k) == 0:
                    continue
                t = float(m[k]['TIME'][0])
                rec = out.setdefault(round(t * 1e9), {})
                for name in ('B', 'V', 'N', 'T'):
                    rec.setdefault(name, []).append(np.asarray(m[k][name], dtype=np.float64))
    return {t: {n: np.concatenate(v) for n, v in rec.items()} for t, rec in out.items()}


base = load(sys.argv[1])
print('%-28s' % 'case' + ''.join('%10s' % s for s in ('t[us]', 'dB', 'dV', 'dN', 'dT')))
for d in sys.argv[2:]:
    c = load(d)
    common = sorted(set(base) & set(c))
    for t in common[2::max(1, len(common) // 6)] + ([common[-1]] if common else []):
        row = '%-28s%10.2f' % (os.path.basename(os.path.normpath(d)), t / 1e3)
        for n in ('B', 'V', 'N', 'T'):
            a, b = base[t][n], c[t][n]
            row += '%10.2e' % (np.linalg.norm(b - a) / max(np.linalg.norm(a), 1e-300))
        print(row)
