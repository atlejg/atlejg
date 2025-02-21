import xarray as xr
import random, glob, os

def rename_and_resample_files(pattern, dpi=300):
    fnms = glob.glob(pattern)
    fig = figure()
    random.shuffle(fnms)
    [os.rename(fnm, f'oppgave_{i+1:02d}.jpg') for i, fnm in enumerate(fnms)]

def showit(p, gridx, nx, ny, save_it=False):
    clf()
    gridy = int(gridx*ny/ny)
    ni, nj = max(int(nx/gridx), 1), max(int(ny/gridy), 1)
    p1 = p0.coarsen(x=ni, y=nj, boundary='trim').mean()
    p1 = p1.reindex_like(p0, method='nearest').astype(uint8)
    imshow(p1.values)
    title(basenm)
    xticks([])
    yticks([])
    tight_layout()
    if save_it:
        fnm1 = f'{basenm}_{i+1:02d}.jpg'
        print(f'creating {fnm1}')
        imsave(fnm1, p1.values)
        return
    show()

fnm0 = sys.argv[1]
basenm, ftype = fnm0.split('.')
assert ftype.lower() == 'jpg'

p = imread(fnm0)
nx, ny, rgb = p.shape
coords = {'x':range(nx), 'y':range(ny), 'rgb':range(rgb)}
p0 = xr.DataArray(p, dims=coords.keys(), coords=coords)

gridxs = [20]
gridxs = [30, 40, 60, 80, 100, 120, 150]

if not 'fig' in dir():
    fig = figure(figsize=(12, 12))
figure(fig.number)

for i, gridx in enumerate(gridxs):
    showit(p, gridx, nx, ny)
    if i == len(gridxs) - 1:
        print('Warning: next picture is not pixelated!')
    ans = input('Press enter to continue (q to jump to full picture:  ')
    if 'q' in ans: break

showit(p, 999999, nx, ny)
