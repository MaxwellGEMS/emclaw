
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for total field intensity and refractive index
    plotfigure = plotdata.new_plotfigure(name='I and n', figno=0)
    
    # plot fields intensity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'I'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = ehfields 
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':2,'markersize':1}

    # plot refractive index

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = '$n(x,t)$'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = refind
    plotitem.plotstyle = '-'
    plotitem.color = 'g'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':2,'markersize':1}

    
    # Figure for q[0], q[1] and n (refractive index)
    plotfigure = plotdata.new_plotfigure(name='E field', figno=1)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = '$E$ and $n$'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = efield
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':1,'markersize':1}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = ''
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = refind
    plotitem.plotstyle = '-'
    plotitem.color = 'g'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':1,'markersize':1}
    
    plotfigure = plotdata.new_plotfigure(name='H field', figno=2)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = '$H$ and $n$'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = hfield
    plotitem.plotstyle = '-'
    plotitem.color = 'r'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':1,'markersize':1}

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = ''
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = refind
    plotitem.plotstyle = '-'
    plotitem.color = 'g'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':1,'markersize':1}

    plotfigure = plotdata.new_plotfigure(name='u and n', figno=3)
    
    # plot fields intensity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'u'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = energy 
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':2,'markersize':1}

    # plot refractive index

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = '$n(x,t)$'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = refind
    plotitem.plotstyle = '-'
    plotitem.color = 'g'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':2,'markersize':1}

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 3           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

def energy(current_data):
    u = 0.5*((current_data.q[0,:]**2)*current_data.aux[0,:] + (current_data.q[1,:]**2)*current_data.aux[1,:])
    return u

def refind(current_data):
    "Return refractive index"
    n = np.sqrt(current_data.aux[0,:]*current_data.aux[1,:])
    return n

def efield(current_data):
    ef = current_data.q[0,:]
    return ef

def hfield(current_data):
    hf = current_data.q[1,:]
    return hf

def ehfields(current_data):
    ehfi = np.sqrt(current_data.q[0,:]**2 + current_data.q[1,:]**2)
    return ehfi

