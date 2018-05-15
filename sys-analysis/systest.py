import sys
sys.path.append('/global/homes/h/huikong/eboss/LSSanalysis/')
sys.path.append('../')
from obiwan_mask import *
from xitools_eboss import *
import healpy as hp
import numpy as np
import astropy.io.fits as fits
from math import sqrt
from optimize import fmin
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages

def Obiwansys(surveyname,map_name='ext',sysmax=0.13,sysmin=0.01,res=256,nest = False,dataset = 'obiwan'):
    filename = surveyname.data
    if dataset == 'obiwan':
    	ranfile = surveyname.obiwan
    elif dataset == 'uniform':
	ranfile = surveyname.uniform
    if map_name == 'star':
	mapname = surveyname.star_map
    elif map_name == 'ext':
	mapname = surveyname.ext_map
    else:
	mapname = surveyname.anand_map
    outputname = './data/'+map_name+'_sys_'+dataset+'.dat'
    print "test1"
    #file readin
    data = fits.open(filename)[1].data
    ran_data = fits.open(ranfile)[1].data
    if map_name == 'ext':
    	fsys = np.loadtxt(mapname)/2.751 #E(B-V)#*2.285/2.751 #r-band extinction
    elif map_name == 'star':
	fsys = np.loadtxt(mapname)
    else:
	fsys = fits.open(mapname)[1].data[map_name]
    #resolution
    npo = 12*res**2
    #galaxies
    pixlg = np.zeros(npo)  
    ng=0.
    print "test2"
    for i in range(0,len(data)):
        ra,dec =  data['ra'][i],data['dec'][i]
        p = int(hp.ang2pix(res,ra,dec,nest = nest,lonlat = True))
        pixlg[p] += 1.
        ng+=1.
    print "test3"
    #randoms
    pixlr = np.zeros(npo)
    nr=0.
    for i in range(0,len(ran_data)):
        ra,dec =  ran_data['ra'][i],ran_data['dec'][i]
        p = int(hp.ang2pix(res,ra,dec,nest = nest,lonlat = True))
        pixlr[p] += 1.#ran_data['obiwan_weight'][i]#modified
        nr += 1.#ran_data['obiwan_weight'][i]#modified
        
    print('total number of galaxies,randoms:')
    print('%f %f' % (ng,nr))
    
    #total bins
    nsysbin = 10
    bing = np.zeros(nsysbin)
    binr = np.zeros(nsysbin)
    
    bintg = 0#total galaxies in good pixel
    bintr = 0#total random in good pixel
    obing = 0#total galaxies outside bin
    obinr = 0#total randoms outside bin
    ng0 = 0#total galaxies in bad pixel
    nr0 = 0#total random in bad pixel
    sysm = float(nsysbin)/(sysmax-sysmin)
    
    for i in range(0,npo):
        sysv = float(fsys[i])
        if sysv != 0: #the maps are not perfect, entries with 0s shouldn't be used
            bintg += pixlg[i]
            bintr += pixlr[i]
            bins = int((sysv-sysmin)*sysm)
            if bins >= 0 and bins < nsysbin:
                bing[bins] += pixlg[i]
                binr[bins] += pixlr[i]
            else:
                obing += pixlg[i]
                obinr += pixlr[i]
        else:
            nr0 += pixlr[i] #count numbers inside bad pixels in sys map
            ng0 += pixlg[i]
            
    print 'total number of randoms/objects '+str(bintr)+'/'+str(bintg)
    print 'number of randoms/objects where sys = 0 '+str(nr0)+'/'+str(ng0)
    print 'number of randoms/objects outside tested range '+str(obinr)+'/'+str(obing)
    ave = float(bintg)/float(bintr)
    print 'average number of objects per random is '+ str(ave)
    
    
    fs = open(outputname,'w')
    
    for i in range(0,nsysbin):
        sysv = sysmin + 1./(2.*sysm) + i/sysm
        if binr[i] > 0:
            ns = bing[i]/binr[i]/ave
            nse = sqrt(bing[i]/(binr[i])**2./(ave)**2.+(bing[i]/ave)**2./(binr[i])**3.)
        else:
            ns = 1.
            nse = 1.
        print('%f %f\n' % (bing[i],binr[i]))
        fs.write(str(sysv)+' '+str(ns)+' '+str(nse)+'\n')
    fs.close()
    print "test4"
    plot_exct(outputname,"./plots/"+map_name+"_"+dataset+".pdf",'eBoss chunk23 '+dataset,xlab = map_name)
    return True

def plot_exct(inputname,output_plot,titlename,xlab = 'ext'):
    print "test5"
    pp = PdfPages(output_plot)
    plt.clf()
    plt.minorticks_on()
    d=np.loadtxt(inputname).transpose()
    chin = sum((d[1]-1.)**2./d[2]**2.)
    print chin
    lf = linfit(d[0],d[1],d[2])
    inl = np.array([1.,0])
    b,m = fmin(lf.chilin,inl)
    print 'b='+str(b)+' m='+str(m)
    chilin = sum((d[1]-(m*d[0]+b))**2./d[2]**2.)
    print chilin
    plt.errorbar(d[0],d[1],d[2],fmt='ko')
    ol = np.ones((len(d[0])))
    plt.plot(d[0],ol,'k:')
    plt.plot(d[0],m*d[0]+b,'k--')
    if xlab == '':
        plt.xlabel(sys,size=16)
    else:
        plt.xlabel(xlab,size=16)

    plt.ylabel(r'$N_{\rm gal}/N_{\rm ran}$ (normalized)',size=16)
    plt.ylim(.7,1.19)
    plt.text(min(d[0])+0.1*(max(d[0])-min(d[0])),1.1,r'$\chi^2$ null ='+str(chin)[:4],color='k')
    plt.text(min(d[0])+0.1*(max(d[0])-min(d[0])),1.08,r'$\chi^2$ lin ='+str(chilin)[:4],color='k')
    #plt.title(r'galaxy density vs. $i$-band depth for v0.7 eboss QSOs, 0.9 < z < 2.2')
    plt.title(titlename)
    pp.savefig()
    pp.close()
    return True   

eboss=surveynames()
eboss.eboss_step5()
#Obiwansys(eboss,'ext',sysmax=0.07,sysmin=0.01)
#Obiwansys(eboss,'star',sysmax=180,sysmin=50)
#Obiwansys(eboss,'hppsfdepth_g',sysmax=1250,sysmin=105,dataset = 'obiwan')
#Obiwansys(eboss,'hppsfdepth_g',sysmax=1250,sysmin=105,dataset = 'uniform')
#Obiwansys(eboss,'hppsfsize_g',sysmax=2.3,sysmin=1.0,dataset = 'obiwan')
#Obiwansys(eboss,'hppsfsize_r',sysmax=2,sysmin=0.9,dataset = 'uniform')
#Obiwansys(eboss,'hppsfsize_r',sysmax=2,sysmin=0.9,dataset = 'obiwan')
#Obiwansys(eboss,'hpstardens',sysmax=2400,sysmin=900,dataset = 'uniform')
#Obiwansys(eboss,'hpstardens',sysmax=2400,sysmin=900,dataset = 'obiwan')
Obiwansys(eboss,'hpebv',sysmax=0.085,sysmin=0.01,dataset = 'uniform')
Obiwansys(eboss,'hpebv',sysmax=0.085,sysmin=0.01,dataset = 'obiwan')

