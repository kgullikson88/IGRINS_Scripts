import numpy
import HelperFunctions
import FittingUtilities
from astropy.io import fits
import astropy.units as u
import sys
import matplotlib.pyplot as plt


def FitWavelength(data, model, fitorder=3):
  data.y /= data.cont
  xgrid = numpy.arange(data.size())
  data_lines = FittingUtilities.FindLines(data, tol=0.98, linespacing=0.2e-3)
  data_flux = data.y[data_lines]
  #data_flux = (data_flux - max(data_flux))/(max(data_flux) - min(data_flux)) + 1

  model.y /= model.cont
  model_lines = FittingUtilities.FindLines(model, tol=0.98, linespacing=0.2e-3)
  model_flux = model.y[model_lines]
  model_lines = model.x[model_lines]


  fig = plt.figure(2, figsize=(10, 7))
  top = fig.add_subplot(211)
  top.plot(data.x, data.y, 'k-')
  top.plot(data.x[data_lines], data_flux, 'ro')
  top.grid(True)
  top.set_ylim((max(-0.1, numpy.min(data.y)), 1.1))
  #top.set_xbound(lower=leftlim, upper=rightlim)

  bottom = fig.add_subplot(212, sharex=top)
  bottom.plot(model.x, model.y, 'k-')
  bottom.plot(model_lines, model_flux, 'ro')
  bottom.grid(True)

  EH = EventHandler(top, bottom, data_lines, data.x[data_lines], model_lines)
  fig.canvas.mpl_connect('button_press_event', EH.OnClick)

  plt.show()

  print EH.pixels
  print EH.waves

  if len(EH.pixels) == len(EH.waves) and len(EH.pixels) >= fitorder:
    fit = numpy.poly1d(numpy.polyfit(EH.pixels, EH.waves, fitorder))
    xaxis = fit(xgrid)
  else:
    xaxis = data.x
  data.y *= data.cont
  return xaxis, EH.pixels, EH.waves




class EventHandler():
  def __init__(self, data_ax, model_ax, data_lines, data_waves, model_lines):
    self.data_ax = data_ax
    self.model_ax = model_ax
    self.pixels = []
    self.waves = []
    self.data_lines = numpy.array(data_lines)
    self.data_waves = numpy.array(data_waves)
    self.model_lines = numpy.array(model_lines)
  

  def OnClick(self, event):
    event.inaxes.plot(event.xdata, event.ydata, 'gx', markersize=10)
    if event.inaxes == self.data_ax:
      #Find the closest data line to this event's xdata
      idx = numpy.argmin(numpy.abs(self.data_waves  - event.xdata))
      self.pixels.append(self.data_lines[idx])
      print "Using line at pixel %i" %self.data_lines[idx]
    elif event.inaxes == self.model_ax:
      #Find the closest model line to this events xdata
      idx = numpy.argmin(numpy.abs(self.model_lines  - event.xdata))
      self.waves.append(self.model_lines[idx])
      print "Using line at wavelength %g um" %self.model_lines[idx]



def CCImprove_Wrapper(data, model, tol=1e-3):
  #Make the data evenly-sampled
  xgrid = numpy.linspace(data.x[0], data.x[-1], data.size())
  data = FittingUtilities.RebinData(data, xgrid)
  dx = xgrid[1] - xgrid[0]

  #Sample the model with the same spacing
  xgrid = numpy.arange(model.x[0], model.x[-1], dx)
  model = FittingUtilities.RebinData(model, xgrid)
  shift = FittingUtilities.CCImprove(data, model, tol=tol, be_safe=True)
  return shift



if __name__ == "__main__":
  import os
  import DataStructures
  telfile = "%s/School/Research/Useful_Datafiles/Telluric_NearIR.dat" %(os.environ["HOME"])
  x,y,c,e = numpy.loadtxt(telfile, unpack=True)
  telluric = DataStructures.xypoint(x=x*u.nm.to(u.micron), y=y, cont=c, err=e)
  telluric = FittingUtilities.ReduceResolution2(telluric, 40000.0)

  for fname in sys.argv[1:]:
    orders = HelperFunctions.ReadExtensionFits(fname)
    columns = []
    for i, order in enumerate(orders):
      if order.x[0] < 1.9:
        band = "H"
      else:
        band = "K"
      wavefile = "%s/School/Research/Useful_Datafiles/IGRINS_Datafiles/SDC%s_order%.3i.wavedata" %(os.environ["HOME"], band, i+2)
      left = numpy.searchsorted(telluric.x, order.x[0] - 2e-3)
      right = numpy.searchsorted(telluric.x, order.x[-1] + 2e-3)
      model = telluric[left:right]

      order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=3, fitorder=2)
      shift = CCImprove_Wrapper(order, model, tol=2e-3)
      print shift
      #plt.plot(ccf.x, ccf.y)
      #plt.show()
      order.x += shift

      done = False
      pixels = []
      waves = []    
      while not done:
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        done = True
        ax.cla()
        ax.plot(order.x, order.y/order.cont, 'k-')
        ax.plot(model.x, model.y, 'r-')
        plt.show(block=False)

        inp = raw_input("Wavelength fit okay? [Y]/n ")
        if "n" in inp.lower():
          done = False
          order.x, newpixels, newwaves = FitWavelength(order, model)

          if len(waves) == 0:
            pixels = newpixels
            waves = newwaves
          else:
            for p, w in zip(newpixels, newwaves):
              if min(abs(numpy.array(pixels) - p)) > 0.1:
                pixels.append(p)
                waves.append(w)

      column = {"wavelength": order.x,
                "flux": order.y,
                "continuum": order.cont,
                "error": order.err}
      columns.append(column)

      if len(pixels) > 3:
        numpy.savetxt(wavefile, numpy.transpose((pixels, waves)))

    HelperFunctions.OutputFitsFileExtensions(column_list, fname, fname, mode="new")




