import matplotlib as mpl
#mpl.use('agg')
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler('color', ['b','g','r','c','m','y','gray','navy','lawngreen','samlon'])
import matplotlib.pyplot as plt
import numpy as np
import gk1 # please update gk1.py in gkeyllscript to most recent version
import os.path

"""
USAGE:
    Copy this file to the data folder, modify following parameter
    as needed, and call "python ohm5m2DAlongY.py" to view/save plots.

    Set save_image to True and pfx to some tag to save the image.
    Set show to True to show image on screen.
"""

Enorm = 1.
iEnorm = 1./Enorm
tnorm = 1.

frames = [6]
filename_base = 'gem'
params = dict(me=1./25., qe=-1., mp=1., qp=1., mo=16., qo=1.,
        num_fluids=3,fluid_names=['e','p','o'],num_moments=5)

save_image = False
pfx = ''
show = True

def calcDtTerm(f, f_fwd, fld, **kwargs):
  idt = 1./(f_fwd.Time() - f.Time())
  v = f.getField(fld, **kwargs)
  v_fwd = f_fwd.getField(fld, **kwargs)
  return (v_fwd - v) * idt

def calcDxTermAlongY(f, fld, ix):
  v = f.getField(fld, ix=slice(ix-1,ix+1))
  nx = f['StructGridField'].shape[0]
  xlo = f.LowerBounds()[0]
  xup = f.UpperBounds()[0]
  dx = (xup-xlo)/nx
  idx = 1./dx
  return (v[:,1]-v[:,0]) * idx

def calcDyTermAlongY(f, fld, ix):
  v = f.getField(fld, ix=ix)
  ret = np.zeros_like(v)
  ny = f['StructGridField'].shape[1]
  ylo = f.LowerBounds()[1]
  yup = f.UpperBounds()[1]
  dy = (yup-ylo)/ny
  idy = 1./dy
  ret[1:] = (ret[1:] - ret[:-1]) * idy
  return ret

def plotOhmAlongY(ts, filename_base, ix=None, species=['e'], hasDtTerm=False):
  f = gk1.File('%s_q_%d.h5'%(filename_base, ts), **params)
  ff = None
  if hasDtTerm:
    ff = gk1.File('%s_f_%d.h5'%(filename_base, ts), **params)
  nx,ny,ncomp = f['StructGridField'].shape
  xlo,ylo = f.LowerBounds()
  xup,yup = f.UpperBounds()
  dx = (xup-xlo)/nx
  dy = (yup-ylo)/ny
  idx = 1./dx
  idy = 1./dy
  if ix == None:
    ix = nx/2
  x = f.getCoordinates('x')
  y = f.getCoordinates('y')
  x_cut = x[ix]
  fig, ax = plt.subplots()
  #
  Ez = f.getField('Ez', ix=ix)
  ax.plot(y, iEnorm*Ez, label='$ E_z $',color='k',lw=3,alpha=0.8)
  #
  for s in species:
    # -VxB term, i.e., convection electric field
    _v_sxB = f.getField('-vx_%s*By+vy_%s*Bx'%(s,s), ix=ix)
    Etotal_s = _v_sxB
    line = ax.plot(y, iEnorm*_v_sxB,
            label='$ -(\\mathbf{v}_%s \\times \\mathbf{B} )|_z$'%(s),
            lw=3,alpha=0.8)[0]
    this_color = line.get_color()
    # div(rho*v*v) term, works for uniform grid only
    inq_s = f.getField('(m%s/q%s)/rho_%s'%(s,s,s), ix=ix)
    div_rhovvz_s = (calcDxTermAlongY(f, 'rhovx_%s*vz_%s'%(s,s), ix=ix)
                  + calcDyTermAlongY(f, 'rhovy_%s*vz_%s'%(s,s), ix=ix)) * inq_s
    ax.plot(y, iEnorm*div_rhovvz_s,
            label='$ \\nabla\\cdot(\\rho_%s\\mathbf{v}_%sv_{z,%s})/{n_%sq_%s}$'%(s,s,s,s,s),
            lw=3,alpha=0.8,color=this_color,ls='--')
    Etotal_s = Etotal_s + div_rhovvz_s
    # d(rho*v)/dt/nq term
    if hasDtTerm:
      drhovz_dt_nq_s =calcDtTerm(f, ff, 'rhovz_'+s, ix=ix) * inq_s
      ax.plot(y, iEnorm*drhovz_dt_nq_s,
              label='$ \\partial_t(\\rho_%s v_{z,%s})/{n_%sq_%s}$'%(s,s,s,s),
              lw=3,alpha=0.8,color=this_color,ls='-.')
      Etotal_s = Etotal_s + drhovz_dt_nq_s
    # sum of different components
    ax.plot(y, iEnorm*Etotal_s,
            label='$ E_{z,%s} $'%(s),
            lw=3,alpha=0.8,color=this_color,ls=':')
  ax.axhline(0, color='gray', lw=2)
  title = '$t=%g, x=%g$'%(f.Time()/tnorm,x_cut)
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title = title, fontsize='x-large')
  lgd.get_title().set_fontsize('x-large')
  ax.set_xlabel('$y$', fontsize='x-large')
  ax.set_ylabel('$E$', fontsize='x-large')
  if save_image:
      fig.savefig('%s_ohm_Ez_ycut_%02d_x%g.png'%(pfx,ts,x_cut), bbox_inches='tight')
  if show:
      plt.show()
  plt.close(fig)
  f.close()
  if ff != None:
      ff.close()

for frame in frames:
  print("frame %d"%frame)
  plotOhmAlongY(frame, filename_base, hasDtTerm=True)
