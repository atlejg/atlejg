import sys, os
import pylab as pl
import math


def create_gridfile(filenm, ny, dy, dx, dz, tops):
    """
    creates a file that defines an adaptive grid for eclipse. by adaptive, i mean that cell size is growing
    outwards from a predfined face (in the xz-plane) that is to contain the well
    (i.e., well is in the y-direction)
    not fully implemented. dx and dz must provided as vectors.
    see create_gridfile2 for a better implementation
    """
    print("creating file %s" % filenm)
    file = open(filenm, 'w')
    #
    file.write('DX\n')
    for k in range(0, len(dz)):
        for j in range(0, ny):
            for i in range(0, len(dx)):
                file.write('%f\n' % dx[i])
    file.write('/\n\n')
    #
    file.write('DY\n')
    file.write('%i*%f\n' % (ny * len(dx) * len(dz), dy))
    file.write('/\n\n')
    #
    file.write('DZ\n')
    for k in range(0, len(dz)):
        file.write('%i*%f\n' % (ny*len(dx), dz[k]))
    file.write('/\n\n')
    file.write('TOPS\n')
    file.write('%i*%f\n' % (ny*len(dx), tops))
    file.write('/\n')
    #
    file.close()
    print('total grid points:', len(dx)*ny*len(dz))


def create_gridfile2(Lx, Ly, Lz, xw, zw, top, dx0, kx, dx_max, dz0, kz, dz_max, ny):
    """
    creates an adaptive grid for eclipse. by adaptive, i mean that cell size is growing
    outwards from a predfined face (in the xz-plane) that is to contain the well.
    #
    note1: well is in the y-direction
    note2: z-direction is downwards
    #
    Lx, Ly, Lz     : width, length, height of reservoir (box)
    xw, zw         : well location in xz-plane
    top            : top of reservoir.
    dx0, dz0       : dimension of center cell
    kx0, kz0       : growth factor outwards from center
    dx_max, dz_max : max cell size
    ny             : number of cells along well (uniform)
    #
    returns a tuple containg:
    a string defining the grid
    the grid dimensions nx, ny, nz
    grid size along the well (dy)
    the indices for the well center (in eclipse coordinates)
    #
    note: you need to set DIMENS manually (it cannot be in GRID-section)
    """
    dy = Ly/float(ny)
    txt, nx,ny,nz, ix,iy,iz = create_gridfile3(Lx, Ly, Lz, xw, 0.5*dy, zw, top, dx0,kx,dx_max, dy,1.,dy, dz0,kz,dz_max)
    return txt, nx,ny,nz, dy, ix,iz

def create_gridfile3(Lx, Ly, Lz, xw, yw, zw, top, dx0, kx, dx_max, dy0, ky, dy_max, dz0, kz, dz_max):
    """
    creates an adaptive grid for eclipse. by adaptive, i mean that cell size is growing
    outwards from a given point
    #
    note: z-direction is downwards
    #
    Li             : width, length, height of reservoir (box). i = x,y,z
    xi             : well location. i = x,y,z
    top            : top of reservoir.
    di0            : dimension of center cell. i = x,y,z
    ki0            : growth factor outwards from center. i = x,y,z
    di_max         : max cell size. i = x,y,z
    #
    returns a tuple containg:
    a string defining the grid
    the grid dimensions nx, ny, nz
    grid size along the well (dy)
    the indices for the well center (in eclipse coordinates)
    #
    note: you need to set DIMENS manually (it cannot be in GRID-section)
    """
    x, ix = get_coord(Lx, xw, dx0, kx, dx_max)
    y, iy = get_coord(Ly, yw, dy0, ky, dy_max)
    z, iz = get_coord(Lz, zw, dz0, kz, dz_max)
    dx, dy, dz = pl.diff(x), pl.diff(y), pl.diff(z)
    nx, ny, nz = len(dx), len(dy), len(dz)
    #
    # build grid-string
    txt  = 'DXV\n'
    for dx_ in dx:
        txt += " %.2f\n" % dx_
    txt += '/\n\n'
    txt += 'DYV\n'
    for dy_ in dy:
        txt += " %.2f\n" % dy_
    txt += '/\n\n'
    txt += 'DZ\n'
    for dz_ in dz:
        txt += " %i*%.2f\n" % (nx*ny, dz_)
    txt += '/\n\n'
    txt += 'BOX\n 1 %i 1 %i 1 1 /\n' % (nx, ny)
    txt += 'TOPS\n %i*%.2f /\n' % (nx*ny, top)
    txt += 'ENDBOX\n'
    #
    return (txt, nx,ny,nz, ix+1,iy+1,iz+1)

def get_coord(Lx, x0, dx1, k, dx_max):
    '''
    make a distribution of points based on a geometric series.
    Lx  = length of line to be gridded
    x0  = where the gridding is finest (typically where the well should be placed)
    dx1 = length of first interval (left and right)
    k   = grow-factor
    dx_max = max length of intervals
    returns an array of grid points and an index for the center.
    '''
    coord = []
    x0 -= dx1/2.      # need this..
    coord.append(x0)
    #
    # going right
    x = x0
    dx = dx1
    while True:
        if dx > dx_max:
            # make sure we dont get the next to last coordinate squeezed into the last.
            # we do this by making homogenous grid on the remaining part
            dL = Lx - x
            n = int(pl.ceil(dL / dx_max))
            dx = dL / float(n)
            for i in range(1, n+1):
                coord.append(x + i*dx)
            break
        x  += dx
        if x >= Lx:
            coord.append(Lx)
            break
        # this will normally happend ...
        coord.append(x)
        dx = dx*k
    #
    # going left
    x = x0
    dx = dx1
    while True:
        dx = dx*k
        if dx > dx_max:
            # make sure we dont get the next to last coordinate squeezed into the last.
            # we do this by making homogenous grid on the remaining part
            dL = x
            n = int(pl.ceil(dL / dx_max))
            dx = dL / float(n)
            for i in range(1, n+1):
                coord.append(x - i*dx)
            break
        x  -= dx
        if x <= 0.:
            coord.append(0.)
            break
        # this will normally happend ...
        coord.append(x)
    #
    coord.sort()
    coord = pl.unique(pl.array(coord))
    center = (coord == x0).nonzero()[0][0]
    return (coord, center)

def tilted_boxgrid(xmax, ymax, zmin, zmax, nx, ny, nz, angle, wedge_factor,
                   directory='.', fname=None, overwrite=False):
    '''
    creates a rectangular box that is tilted downwards in the x-direction.
    it uses COORD to define the xy-coordinates and ZCORN to define the depths.
    this one does not use BOX
    if no filename is provided, it will create one based on input parameters.
    the file will be created in the given directory.
    [length] = m
    [angle]  = deg. this is the angle of the top reservoir surface.
    wedge_factor: a value of 1 gives no wedging, 0.5 means thickness of the reservoir on the
                  shallow side is half of the thickness on the deep side.
    '''
    # check fname and if gridfile already exists
    if fname is None:
        fname = '%.1fx%.1fx%.1f-%.1f_%ix%ix%i_%.1f_%.2f.GRDECL' % \
              (xmax, ymax, zmin, zmax, nx, ny, nz, angle, wedge_factor)
    fname = '%s/%s' % (directory, fname)
    if os.path.exists(fname) and not overwrite:
        print('%s already exists' % fname)
        return fname
    #
    # so - we will have to create the grid.
    #
    # calculate the height variance from first to last row in x direction
    dhmax = pl.sin(math.radians(angle)) * xmax if angle > 0 else 0.
    dh = dhmax * pl.linspace(0, 1, nx+1)
    # varying dz for wedge
    dz_scaler = pl.linspace(wedge_factor, 1., nx+1)
    # project in the dipping direction
    xmax *= pl.cos(math.radians(angle))
    #
    dy = ymax / float(ny)
    dx = xmax / float(nx)
    dz = (zmax-zmin) / float(nz)
    #
    f = open(fname, 'w')
    #
    f.write('''-- created by tilted_boxgrid()
 -- input parameters:
 -- xmax = %.1f ymax = %.1f zmin = %1.f zmax = %.1f
 -- nx = %i ny = %i nz = %i
 -- angle = %.1f deg
 -- wedge_factor = %.3f
 -- filename = %s''' %
       (xmax, ymax, zmin, zmax, nx, ny, nz, angle, wedge_factor, fname)
    )
    #
    f.write('\n\nCOORD\n')
    for j in pl.arange(ny+1):
        for i in pl.arange(nx+1):
            x = i*dx
            y = j*dy
            f.write('%.1f %.1f %.1f   %.1f %.1f %.1f\n' % (x,y,0, x,y,0))
    f.write('/\n\n')
    #
    f.write('ZCORN\n')
    #
    # indices should be like (0, 1,1, 2,2, 3,3, ... , nx/ny/nz)
    i_indices = list(range(nx)) + list(pl.arange(nx)+1)
    i_indices.sort()
    j_indices = list(range(ny)) + list(pl.arange(ny)+1)
    j_indices.sort()
    k_indices = list(range(nz)) + list(pl.arange(nz)+1)
    k_indices.sort()
    for k in k_indices:
        for j in j_indices:
            for i in i_indices:
                h = zmin + k*dz*dz_scaler[i] + dh[i]
                f.write('%.3f\n' % h)
            f.write('\n')
        f.write('\n--\n')
    f.write('/\n\n')
    f.close()
    print('gridfile %s created' % fname)
    return fname

if __name__ == '__main__':

    '''
    # simulation box
    Lx = 1500
    Ly = 3000
    Lz = 45
    tops = 1738

    # well-center. we specify the xz-face for the cell where the well recides. horiz part of well is in y-dir.
    # given by upper left corner coordinates and the length of sides
    x0  = 749
    z0  = 17.5 # relative to upper limit (tops)
    dx1 = 2
    dz1 = 0.2

    # grid parameters. growth factors and limits
    kx = 1.2
    kz = 1.2
    dx_max = 50
    dz_max = 1.0
    dy = 100

    # calculations
    ny = Ly / dy
    (x_coord, i) = get_coord(Lx, x0, dx1, kx, dx_max)
    (z_coord, k) = get_coord(Lz, z0, dz1, kz, dz_max)
    dx = pl.diff(x_coord)
    dz = pl.diff(z_coord)

    directory = sys.argv[1]

    filenm = 'grid_%ix%ix%i_%ix%ix%i_%.1f_%.1f_%1.f_%.1f_%.1f_%i_%.1f_.inc' % \
       (len(dx), ny, len(dz), i, int(ny/2), k, dx1, dy, dz1, kx, kz, dx_max, dz_max)
    create_gridfile('%s/%s' % (directory, filenm), ny, dy, dx, dz, tops)
    '''
    tilted_boxgrid(1,1,0,1, 3,3,2, 45, overwrite=True)
