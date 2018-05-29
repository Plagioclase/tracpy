"""
Plotting routines for tracking.
"""

from __future__ import absolute_import
# import matplotlib as mpl
# # set matplotlib to use the backend that does not require a windowing system
# mpl.use("Agg")
import numpy as np
# from matplotlib.mlab import *
import matplotlib.pyplot as plt
from . import inout
import os
import matplotlib.ticker as ticker
import matplotlib as mpl
from . import op
from . import tools
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature
import cmocean.cm as cmo

pc = ccrs.PlateCarree()

def background(grid, proj=ccrs.Mercator(), ax=None, pars=np.arange(18, 35),
               mers=np.arange(-100, -80), figsize=None,
               hlevs=np.hstack(([10, 20], np.arange(50, 500, 50))),
               extent=[-98, -87.5, 22.8, 30.5],
               col='lightgrey', halpha=1, fig=None, outline=[1, 1, 0, 1],
               res='50m', plotstates=True, fontsize=14):
    """
    Plot basic TXLA shelf background: coastline, bathymetry, meridians, etc
    Can optionally input grid (so it doesn't have to be loaded again)

    Args:
        proj
        extent: west lon, east lon, south lat, north lat for axes bounds
        pars: parallels to plot
        mers: meridians to plot
        hlevs: which depth contours to plot
        outline: west, east, north, south lines (left, right, top, bottom)
        res: 10m, 50m, 110m
    """

    mpl.rcParams.update({'font.size': fontsize})#,'font.weight': 'bold'})
    land = cfeature.NaturalEarthFeature('physical', 'land', res,
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
    states = cfeature.NaturalEarthFeature('cultural',
                                          name='admin_1_states_provinces_lines',
                                          scale=res, facecolor='none')

    if fig is None and figsize is None:
        fig = plt.figure(figsize=(9.4, 7.7), dpi=100)
    elif fig is None and figsize is not None:
        fig = plt.figure(figsize=figsize, dpi=100)

    if ax is None:
        ax = fig.add_subplot(111, projection=proj)

    ax.set_extent(extent, pc)
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_bottom = False  # turn off labels where you don't want them
    gl.ylabels_right = False

    # plot isobaths
    ax.contour(grid.lon_rho, grid.lat_rho, grid.h, hlevs, colors='0.6',
               transform=pc, linewidths=0.5, alpha=halpha)

    # Do plot
    ax.add_feature(land, facecolor='0.8')
    ax.coastlines(resolution=res)  # coastline resolution options are '110m', '50m', '10m'
    if plotstates:
        ax.add_feature(states, edgecolor='0.2')
    ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='0.2')

    # Outline numerical domain
    # if outline:  # backward compatibility
    #     outline = [1,1,1,1]
    if outline[0]:  # left
        ax.plot(grid.lon_rho[:, 0], grid.lat_rho[:, 0], 'k:', transform=pc)
    if outline[1]:  # right
        ax.plot(grid.lon_rho[:, -1], grid.lat_rho[:, -1], 'k:', transform=pc)
    if outline[2]:  # top
        ax.plot(grid.lon_rho[-1, :], grid.lat_rho[-1, :], 'k:', transform=pc)
    if outline[3]:  # bottom
        ax.plot(grid.lon_rho[0, :], grid.lat_rho[0, :], 'k:', transform=pc)

    return fig, ax


def hist(xp, yp, proj, fname, grid, tind='final', which='contour', vmax=None,
         fig=None, ax=None, bins=(40, 40), N=10, xlims=None,
         ylims=None, C=None, Title=None, weights=None, cmap=None,
         Label='Final drifter location (%)', binscale=None,
         cbcoords=[0.35, 0.25, 0.6, 0.02],
         crsproj=ccrs.LambertConformal()):
    """
    update
    Plot histogram of given track data at time index tind.

    Args:
        xp,xp: Drifter track positions in projected coordinates [ndrifters x time]
        fname: Plot name to save
        tind (Optional): Default is 'final', in which case the final
         position of each drifter in the array is found and plotted.
         Alternatively, a time index can be input and drifters at that time
         will be plotted. Note that once drifters hit the outer numerical
         boundary, they are nan'ed out so this may miss some drifters.
        which (Optional[str]): 'contour', 'pcolor', 'hexbin', 'hist2d' for
         type of plot used. Default 'hexbin'.
        bins (Optional): Number of bins used in histogram. Default (15,25).
        N (Optional[int]): Number of contours to make. Default 10.
        grid (Optional): grid as read in by inout.readgrid()
        xlims (Optional): value limits on the x axis
        ylims (Optional): value limits on the y axis

    Note:
        Currently assuming we are plotting the final location of each drifter
        regardless of tind.

    Example usage:
        import tracpy
        import tracpy.plotting
        import netCDF4 as netCDF
        proj = tracpy.tools.make_proj('nwgom-pyproj')
        grid = tracpy.inout.readgrid([gridlocation], proj)
        nc = netCDF.Dataset([trackfile.gc.nz])
        xg = nc['xg'][:]; yg = nc['yg'][:]
        xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
        tracpy.plotting.hist(xp, yp, fname, grid)
    """

    # # Change positions from lon/lat to x/y for analysis
    # xp, yp = grid.proj(lonp, latp)
    # # Need to retain nan's since basemap changes them to values
    # ind = np.isnan(lonp)
    # xp[ind] = np.nan
    # yp[ind] = np.nan

    if fig is None:
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111)

    if tind == 'final':
        # Find final positions of drifters
        xpc, ypc = tools.find_final(xp, yp)
    elif isinstance(tind, int):
        xpc = xp[:, tind]
        ypc = yp[:, tind]
    else:  # just plot what is input if some other string
        xpc = xp.flatten()
        ypc = yp.flatten()

    if cmap is None:
        cmap = 'YlOrRd'

    if which == 'contour':

        # Info for 2d histogram
        H, xedges, yedges = np.histogram2d(xpc, ypc,
                                           range=[[grid.x_rho.min(),
                                                   grid.x_rho.max()],
                                                  [grid.y_rho.min(),
                                                   grid.y_rho.max()]],
                                           bins=bins)

        # Contour Plot
        XE, YE = np.meshgrid(op.resize(xedges, 0), op.resize(yedges, 0))
        d = (H/H.sum())*100
        # # from http://matplotlib.1069221.n5.nabble.com/question-about-contours-and-clim-td21111.html
        # locator = ticker.MaxNLocator(50) # if you want no more than 10 contours
        # locator.create_dummy_axis()
        # locator.set_bounds(0,1)#d.min(),d.max())
        # levs = locator()
        con = ax.contourf(XE, YE, d.T, N, transform=crsproj)  # ,levels=levs)#(0,15,30,45,60,75,90,105,120))
        con.set_cmap(cmap)

        if Title is not None:
            ax.set_title(Title)

        # Horizontal colorbar below plot
        cax = fig.add_axes(cbcoords)  # colorbar axes
        cb = fig.colorbar(con, cax=cax, orientation='horizontal')
        cb.set_label('Final drifter location (percent)')

        # Save figure into a local directory called figures. Make directory
        # if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        fig.savefig('figures/' + fname + 'histcon.png', bbox_inches='tight')

    elif which == 'pcolor':

        # Info for 2d histogram
        H, xedges, yedges = np.histogram2d(xpc, ypc,
                                           range=[[grid.x_rho.min(),
                                                   grid.x_rho.max()],
                                                  [grid.y_rho.min(),
                                                   grid.y_rho.max()]],
                                           bins=bins, weights=weights)

        # Pcolor plot
        # C is the z value plotted, and is normalized by the total number of
        # drifters
        if C is None:
            C = (H.T/H.sum())*100
        else:
            # or, provide some other weighting
            C = (H.T/C)*100

        p = ax.pcolormesh(xedges, yedges, C, cmap=cmap, transform=crsproj)

        if Title is not None:
            ax.set_title(Title)

        # Set x and y limits
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)

        # Horizontal colorbar below plot
        cax = fig.add_axes(cbcoords)  # colorbar axes
        cb = fig.colorbar(p, cax=cax, orientation='horizontal')
        cb.set_label('Final drifter location (percent)')

        # Save figure into a local directory called figures. Make directory
        # if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        fig.savefig('figures/' + fname + 'histpcolor.png', bbox_inches='tight')
        # savefig('figures/' + fname + 'histpcolor.pdf',bbox_inches='tight')

    elif which == 'hexbin':

        if ax is None:
            ax = plt.gca()
        else:
            ax = ax

        if C is None:
            # C with the reduce_C_function as sum is what makes it a percent
            C = np.ones(len(xpc))*(1./len(xpc))*100
        else:
            C = C*np.ones(len(xpc))*100
        hb = ax.hexbin(xpc, ypc, C=C, cmap=cmap, gridsize=bins[0],
                    extent=(grid.x_psi.min(), grid.x_psi.max(),
                            grid.y_psi.min(), grid.y_psi.max()),
                    reduce_C_function=sum, vmax=vmax, axes=ax, bins=binscale,
                    transform=crsproj)

        # Set x and y limits
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)

        if Title is not None:
            ax.set_title(Title)

        # Want colorbar at the given location relative to axis so this works
        # regardless of # of subplots, so convert from axis to figure
        # coordinates. To do this, first convert from axis to display coords
        # transformations:
        # http://matplotlib.org/users/transforms_tutorial.html
        # axis: [x_left, y_bottom, width, height]
        # ax_coords = [0.35, 0.25, 0.6, 0.02]
        # display: [x_left,y_bottom,x_right,y_top]
        disp_coords = ax.transAxes.transform([(cbcoords[0], cbcoords[1]),
                                              (cbcoords[0]+cbcoords[2],
                                               cbcoords[1]+cbcoords[3])])
        # inverter object to go from display coords to figure coords
        inv = fig.transFigure.inverted()
        # figure: [x_left,y_bottom,x_right,y_top]
        fig_coords = inv.transform(disp_coords)
        # actual desired figure coords. figure:
        # [x_left, y_bottom, width, height]
        fig_coords = [fig_coords[0, 0], fig_coords[0, 1], fig_coords[1, 0] -
                      fig_coords[0, 0], fig_coords[1, 1] - fig_coords[0, 1]]
        # Inlaid colorbar
        cax = fig.add_axes(fig_coords)

        # # Horizontal colorbar below plot
        # cax = fig.add_axes([0.3775, 0.25, 0.48, 0.02]) # colorbar axes
        cb = fig.colorbar(hb, cax=cax, orientation='horizontal')
        cb.set_label(Label)

        # Save figure into a local directory called figures. Make directory
        # if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        fig.savefig('figures/' + fname + 'histhexbin.png', bbox_inches='tight')
        # savefig('figures/' + fname + 'histhexbin.pdf',bbox_inches='tight')

    elif which == 'hist2d':

        ax.hist2d(xpc, ypc, bins=40, range=[[grid.x_rho.min(),
                                          grid.x_rho.max()],
                                         [grid.y_rho.min(),
                                          grid.y_rho.max()]], normed=True)
        plt.set_cmap(cmap)
        # Set x and y limits
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)

        # Horizontal colorbar below plot
        cax = fig.add_axes(cbcoords)  # colorbar axes
        cb = fig.colorbar(cax=cax, orientation='horizontal')
        cb.set_label('Final drifter location (percent)')

        # Save figure into a local directory called figures. Make directory
        # if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        fig.savefig('figures/' + fname + 'hist2d.png', bbox_inches='tight')
        # savefig('figures/' + fname + 'histpcolor.pdf',bbox_inches='tight')


def tracks(lonp, latp, fname, grid, fig=None, ax=None, Title=None):
    """
    Plot tracks as lines with starting points in green and ending points in
    red.

    Args:
        lonp,latp: Drifter track positions [time x ndrifters]
        fname: Plot name to save
    """

    if fig is None:
        fig = plt.figure(figsize=(12, 10))

    if ax is None:
        ax = plt.gca()

    # Starting marker
    ax.plot(lonp[:, 0], latp[:, 0], 'o', color='g', markersize=3,
            label='_nolegend_', alpha=0.4, transform=pc)

    # Plot tracks
    ax.plot(lonp.T, latp.T, '-', color='grey', linewidth=.2, transform=pc)

    # Find final positions of drifters
    lonpc, latpc = tools.find_final(lonp, latp)
    ax.plot(lonpc, latpc, 'o', color='r', label='_nolegend_', transform=pc)

    if Title is not None:
        ax.set_title(Title)

    # Legend, of sorts
    # ax = gca()
    xtext = 0.45
    ytext = 0.18
    ax.text(xtext, ytext, 'starting location', fontsize=16, color='green',
         alpha=.8, transform=ax.transAxes)
    ax.text(xtext, ytext-.03, 'track', fontsize=16, color='grey',
         transform=ax.transAxes)
    ax.text(xtext, ytext-.03*2, 'ending location', fontsize=16, color='red',
         transform=ax.transAxes)

    # Save figure into a local directory called figures. Make directory if it
    # doesn't exist.
    if not os.path.exists('figures'):
        os.makedirs('figures')

    fig.savefig('figures/' + fname + 'tracks.png', bbox_inches='tight')


def transport(name, fmod=None, Title=None, dmax=None, N=7, extraname=None,
              llcrnrlon=-98.5, llcrnrlat=22.5, urcrnrlat=31.0,
              urcrnrlon=-87.5, colormap='Blues', fig=None, ax=None):
    """
    Make plot of zoomed-in area near DWH spill of transport of drifters over
    time.

    FILL IN

    Args:
        name
        U
        V
        lon0
        lat0
        T0
    """

    # Load in transport information
    U, V, lon0, lat0, T0 = inout.loadtransport(name, fmod=fmod)

    # Smaller basemap parameters.
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = inout.readgrid(loc, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                          urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon)

    S = np.sqrt(op.resize(U, 1)**2+op.resize(V, 0)**2)
    Splot = (S/T0)*100
    if dmax is None:
        dmax = Splot.max()
    else:
        dmax = dmax
    # from http://matplotlib.1069221.n5.nabble.com/question-about-contours-and-clim-td21111.html
    locator = ticker.MaxNLocator(N)  # if you want no more than 10 contours
    locator.create_dummy_axis()
    locator.set_bounds(0, dmax)  # d.min(),d.max())
    levs = locator()

    if fig is None:
        fig = plt.figure(figsize=(11, 10))
    else:
        fig = fig
    background(grid=grid)
    c = fig.contourf(grid.xpsi, grid.ypsi, Splot, cmap=cmap,
                 extend='max', levels=levs)
    plt.title(Title)

    # # Add initial drifter location (all drifters start at the same location)
    # lon0 = lon0.mean()
    # lat0 = lat0.mean()
    # x0, y0 = grid['basemap'](lon0, lat0)
    # plot(x0, y0, 'go', markersize=10)

    if ax is None:
        ax = plt.gca()
    else:
        ax = ax
    # Want colorbar at the given location relative to axis so this works
    # regardless of # of subplots,
    # so convert from axis to figure coordinates
    # To do this, first convert from axis to display coords
    # transformations: http://matplotlib.org/users/transforms_tutorial.html
    ax_coords = [0.35, 0.25, 0.6, 0.02]  # axis: [x_left, y_bottom, width, height]
    # display: [x_left,y_bottom,x_right,y_top]
    disp_coords = ax.transAxes.transform([(ax_coords[0], ax_coords[1]),
                                          (ax_coords[0]+ax_coords[2],
                                           ax_coords[1]+ax_coords[3])])
    # inverter object to go from display coords to figure coords
    inv = fig.transFigure.inverted()
    # figure: [x_left,y_bottom,x_right,y_top]
    fig_coords = inv.transform(disp_coords)
    # actual desired figure coords. figure: [x_left, y_bottom, width, height]
    fig_coords = [fig_coords[0, 0], fig_coords[0, 1], fig_coords[1, 0] -
                  fig_coords[0, 0], fig_coords[1, 1] - fig_coords[0, 1]]

    # Inlaid colorbar
    cax = fig.add_axes(fig_coords)
    # cax = fig.add_axes([0.39, 0.25, 0.49, 0.02])
    # cax = fig.add_axes([0.49, 0.25, 0.39, 0.02])
    cb = fig.colorbar(cax=cax, orientation='horizontal')
    cb.set_label('Normalized drifter transport (%)')

    if extraname is None:
        fig.savefig('figures/' + name + '/transport', bbox_inches='tight')
    else:
        fig.savefig('figures/' + name + '/' + extraname + 'transport',
                bbox_inches='tight')
