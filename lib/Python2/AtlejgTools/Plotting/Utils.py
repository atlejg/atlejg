import pylab as pl
import pdb

def plot_scattered_data(x, y, z, nx=100, ny=100, plot_dots=True):
   '''
   plots 2d data.
   x, y : coordinate
   z    : data value
   based on scipy cookbook (http://www.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data)
   '''
   import numpy as np
   from matplotlib.mlab import griddata
   import matplotlib.pyplot as plt
   import numpy.ma as ma
   from numpy.random import uniform
   # make up some randomly distributed data
   # define grid.
   xi = np.linspace(x.min(), x.max(), nx)
   yi = np.linspace(y.min(), y.max(), ny)
   # grid the data.
   zi = griddata(x,y,z,xi,yi)
   #
   # contour the gridded data, plotting dots at the randomly spaced data points.
   plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
   plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
   plt.colorbar() # draw colorbar
   # plot data points.
   if plot_dots: plt.scatter(x,y,marker='o',c='b',s=5)
   #plt.xlim(x.min(), x.max())
   #plt.ylim(y.min(), y.max())
   plt.show()

def contourf(x, y, z, zmin=None, zmax=None, fig=None, nlvls=20, **kwargs):
   '''
   Make contourplot with a useful colormap scaled to zmin, zmax
   Based on http://bytes.com/topic/python/answers/891999-matplotlib-colorbar-scale-problem
   Ex:
      ax = PU.contourf(t, x, z/maxval, 0., 1.)
      setp(ax, 'xlabel', 'TIME')
      setp(ax, 'ylabel', 'MD')
      setp(ax, 'title', '%s %s (relative)' % (frftfile.split('.')[0], varnm))
   '''
   if fig is None: fig = pl.figure()
   if zmin is None: zmin = z.ravel().min()
   if zmax is None: zmax = z.ravel().max()
   fig.clf()
   norm = pl.Normalize(vmin=zmin, vmax=zmax)
   #
   # plot contours
   ax1 = fig.add_axes([0.11, 0.09, 0.73, 0.81], kwargs) # kwargs does not seem to work
   pl.plt.contourf(x, y, z, nlvls, norm=norm)
   pl.plt.contourf(x, y, z, nlvls, norm=norm) # plotting twice to avoid contourlines. bug??
   pl.axis((x[0], x[-1], y[0], y[-1]))    # not sure why this is necessary
   #
   # plot colorbar
   ax2 = fig.add_axes([0.86, 0.09, 0.06, 0.81])
   pl.mpl.colorbar.ColorbarBase(ax2, norm=norm, orientation='vertical')
   #
   return ax1

def contours(z, zmin=None, zmax=None, fig=None, nlvls=20):
   x = pl.arange(z.shape[0])
   y = pl.arange(z.shape[1])
   return contourf(x, y, z.T, zmin=zmin, zmax=zmax, fig=fig, nlvls=nlvls)
