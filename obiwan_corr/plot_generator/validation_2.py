eboss_2chunks = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/corrfunc_data/eboss_final041918__5degree.out'
Johan_random = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/corrfunc_data/Johan_random.txt'
uniform_random = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/corrfunc_data/uniform_eboss_final041918__5degree.out'
weighed_uniform = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/corrfunc_data/weighed_uniform.out'
import numpy as np
import matplotlib.pyplot as plt
ran_eboss_dat = np.loadtxt(eboss_2chunks).transpose()
Johan_random = np.loadtxt(Johan_random).transpose()
uniform_random = np.loadtxt(uniform_random).transpose()
weighed_random = np.loadtxt(weighed_uniform).transpose()
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
major_ticks = np.array([0.01,0.1,1])
minor_ticks = np.arange(0, 0.3, 0.01)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)

ylim = 0.026
major_ticks = np.arange(0,ylim,0.002)
minor_ticks = np.arange(0, ylim, 0.001)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)
ax.grid(which='minor', alpha=0.2,linestyle='--')
ax.grid(which='major', alpha=0.8,linestyle='--')

ax.set_xlim((0.009,5))
ax.set_ylim((0,ylim))

obiwan = plt.errorbar(ran_eboss_dat[0],ran_eboss_dat[1]*ran_eboss_dat[0],ran_eboss_dat[2]*ran_eboss_dat[0],label = 'obiwan')
Johan = plt.errorbar(Johan_random[0],Johan_random[1]*Johan_random[0],Johan_random[2]*Johan_random[0],label = 'Johan')
uniform = plt.errorbar(uniform_random[0],uniform_random[1]*uniform_random[0],uniform_random[2]*uniform_random[0],label = 'uniform')
weighed = plt.errorbar(weighed_random[0],weighed_random[1]*weighed_random[0],weighed_random[2]*weighed_random[0],label = 'weighed')
plt.legend()
plt.xlabel(r'$\theta$'+'(Degree)')
plt.ylabel(r'$\theta$'+'*w('+r'$\theta$'+')')
plt.title('eBoss chunk23')
plt.gca().set_xscale("log", nonposx='clip')
from matplotlib.backends.backend_pdf import PdfPages
with PdfPages('./plots/validation2.pdf') as pdf:
	pdf.savefig() 
plt.show()
