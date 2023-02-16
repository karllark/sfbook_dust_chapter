import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

from astropy.table import Table

# setup the plots
fontsize = 16
font = {'size'   : fontsize}

mpl.rc('font', **font)

mpl.rc('lines', linewidth=2)
mpl.rc('axes', linewidth=2)
mpl.rc('xtick.major', width=2)
mpl.rc('ytick.major', width=2)
mpl.rc('xtick.minor', width=2)
mpl.rc('ytick.minor', width=2)

fig, ax = plt.subplots(figsize=(10.,10.0))
gs = gridspec.GridSpec(2, 2)
ax = []
ax.append(plt.subplot(gs[0,1]))
ax.append(plt.subplot(gs[0,0]))
ax.append(plt.subplot(gs[1,1]))
ax.append(plt.subplot(gs[1,0]))

filetag = ['1','2','3','4','5']
av1 = np.array([0.0, 0.2, 0.5, 1.0, 2.0])

files = ['dirty/DIRTY_bm1_effgrain_tau_001.00_theta_000_mod' 
         + k + '_global_lum.table.fits' for k in filetag]

# get the av=2.0 embedded case
a = Table.read(files[4])
x = 1.0/a['wavelength']
print(a.colnames)
# get the index for the "V" band
sindxs = np.argsort(np.abs(x - (1.0/0.55)))
    
# now mix with the other cases
files = files[0:4]
for k, cfile in enumerate(files):
    b = Table.read(cfile)

#    c_input = b['Flux_Input']
#    c_output_direct = b['Flux_rt_d']
#    c_output_scat = b['Flux_rt_s']
#    c_output = b['Flux']

    c_input = a['Flux_Input'] + b['Flux_Input']
    c_output_direct = a['Flux_rt_d'] + b['Flux_rt_d']
    c_output_scat = a['Flux_rt_s'] + b['Flux_rt_s']
    c_output = a['Flux'] + b['Flux']
    
    tau_eff = -2.5*np.log10(c_output/c_input)
    tau_eff_direct = -2.5*np.log10(c_output_direct/c_input)

    print(c_output_direct[sindxs[0]],c_input[sindxs[0]])
    
    #print(a['Flux_rt_s'], a['Flux_Input'])
    #ax.plot(1./x, a['Flux_Input'], label='mod_i' + str(k))
    #ax.plot(1./x, a['Flux_rt_d'], label='mod_d' + str(k))

    label_str = 'Model{:1d}, $Att(V) = {:3.2f}$'
    ax[0].plot(1./x, tau_eff/tau_eff[sindxs[0]], 
            label=label_str.format(k+1,tau_eff[sindxs[0]]))

    label_str = 'Model{:1d}, $Att(V) = {:3.2f}$'
    ax[1].plot(1./x, tau_eff_direct/tau_eff_direct[sindxs[0]], 
            label=label_str.format(k+1,tau_eff_direct[sindxs[0]]))

    label_str = '$A_1(V) = {:3.1f}$, $Att(V) = {:3.2f}$'
    ax[2].plot(1./x, c_output_scat/c_input, 
            label=label_str.format(av1[k],tau_eff_direct[sindxs[0]]))

ax[0].plot(1./x, a['tau_norm']/a['tau_norm'][sindxs[0]], 
           'k--', label='input extinction')
ax[1].plot(1./x, a['tau_norm']/a['tau_norm'][sindxs[0]], 
           'k--', label='input extinction')

ax[2].set_xlabel('$\lambda$ [$\mu m$]')
#ax[3].set_xlabel('$\lambda$ [$\mu m$]')
ax[0].set_ylabel('$Att(\lambda)/Att(V)$')
ax[1].set_ylabel('$Att(\lambda)/Att(V)$')
ax[2].set_ylabel('$F(scat)/F(total)$')

ax[0].yaxis.tick_right()
ax[0].yaxis.set_label_position("right")

ax[2].yaxis.tick_right()
ax[2].yaxis.set_label_position("right")

#ax.set_yscale('log')

ax[0].legend(loc='best', fontsize=12, title='Extinction + Scattering')
ax[1].legend(loc='best', fontsize=12, title='Extinction Only')

for i in range(3):

    ax[i].tick_params(which='major', length=7)
    ax[i].tick_params(which='minor', length=4)

    ax[i].set_xscale('log')
    ax[i].set_xlim(1./10.,1.0/0.3)

ax[0].set_ylim(0.0,6.5)
ax[1].set_ylim(0.0,6.5)

# now draw the model geometry

ax[3].plot([-2.0,-5.0,-5.0,-2.0,-2.0],
           [-5.0,-5.0,5.0,5.0,-5.0], 'k--')
ann_x = [-2.0,-2.3, -2.75, -3.5]
ann_y = [0.0,0.0,0.0,0.0]
ax[3].plot(ann_x, ann_y, 'ro')
for k in range(len(ann_x)):
    if k % 2 == 0:
        yval = 3.0
        valign = 'center'
    else:
        yval = -3.0
        valign = 'center'
    ax[3].annotate('Star #1, Model' + str(k+1), 
                   xy=(ann_x[k], ann_y[k]), xytext=(ann_x[k]+2.0,yval), 
                   rotation=90., verticalalignment=valign, fontsize=12., 
                   arrowprops=dict(facecolor='black', shrink=0.1,
                                   width=2))

ax[3].plot([-5.0],[0.0], 'bo')
ax[3].annotate('Star #2', xy=(-5., 0.0), xytext=(-4.0,2.0), 
               rotation=90., verticalalignment="center", fontsize=12., 
               arrowprops=dict(facecolor='black', shrink=0.1,
                               width=2))

ax[3].plot([5.0,3.8,4.0,3.85,3.8,3.85,4.0,3.8,5.0],
           [0.0,1.1,0.925,0.5,0.0,-0.5,-0.925,-1.1,0.0],
           'k-')
ax[3].plot([3.85,4.05,4.1,4.05,3.85],
           [0.25,0.125,0.0,-0.125,-0.25],'k-')

ax[3].annotate('Observer', xy=(4.4, 1.0), xytext=(4.0,3.0), 
               rotation=0., verticalalignment="center", 
               horizontalalignment="center",
               fontsize=12., 
               arrowprops=dict(facecolor='black', shrink=0.1,
                               width=2))

ax[3].set_xlim(-5.5,5.5)
ax[3].set_ylim(-5.5,5.5)
ax[3].axis('off')

fig.tight_layout()

fig.savefig('2star_mix_dirty.pdf')
