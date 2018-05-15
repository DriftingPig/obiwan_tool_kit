'''
split the files into multiple piece to run parelle angular correlation function
'''
import sys
sys.path.append('../')
from obiwan_mask import *
import healpy as hp
from astropy.table import Table
from math import *
dirname = "weighed_uniform"
from subprocess import call
call(["mkdir", "/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/splitdata/"+dirname+"/"])
output_filename = "/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/splitdata/"+dirname+"/"
def radec2thphi(ra,dec):
        return (-dec+90.)*pi/180.,ra*pi/180.
    
def ranHelpsort(filename,res=256,rad=''):
    origin_dat = fits.open(filename)[1].data
    npix = 12*res**2
    angm = 1.
    if rad == 'rad':
        angm = 180./pi
    pixls=np.zeros(len(origin_dat))
    print len(origin_dat)
    for i in range(0,len(origin_dat)):
        ra,dec = float(origin_dat['ra'][i])*angm,float(origin_dat['dec'][i])*angm
        th,phi = radec2thphi(ra,dec)
        pixls[i]=hp.pixelfunc.ang2pix(res,th,phi,nest=True)
    print 'done'
    orig_cols = origin_dat.columns
    pixls = np.array(pixls)
    new_cols = fits.ColDefs([
            fits.Column(name='pixl',format='I',
                        array=pixls)
                        ])
    tb_hdu = fits.BinTableHDU.from_columns(orig_cols+new_cols)
    tb_dat = Table.read(tb_hdu)
    tb_dat.sort('pixl')
    return tb_dat

#generate txt files which inculde sin, cos of ra,dec
def SprtFile_nest(filename, total_subs,res=256,rad=''):
    dat = ranHelpsort(filename,res,rad)
    import os
    base_filename = os.path.basename(filename)
    files = []
    for i in range(0,total_subs):
        sub_filename = base_filename[:-5]+'_subset'+str(i)+'.dat'
        f = open(os.path.join(output_filename,sub_filename),'w')
        files.append(f)
    
    sub_num = len(dat)/total_subs
    count = 0
    j=0
    pixflag=[]
    for i in range(0,len(dat)):
        ra = dat['ra'][i]
        dec = dat['dec'][i]
        thi = dec*pi/180.
        phi = ra*pi/180.
        files[j].write(str(sin(thi))+' '+str(cos(thi))+' '+str(sin(phi))+' '+str(cos(phi))+' '+str(dat['obiwan_weight'][i])+'\n')
        if j != total_subs-1:
            count+=1
            if count == sub_num:
                pixflag.append(dat['pixl'][i])
                j+=1
                count=0
    for i in range(0,total_subs):
        files[i].close()
    return pixflag

def SprtFile_nest_dat(filename, total_subs, pixflag, res=256,rad=''):
    dat = ranHelpsort(filename,res,rad)
    import os
    base_filename = os.path.basename(filename)
    files = []
    for i in range(0,total_subs):
        sub_filename = base_filename[:-5]+'_subset'+str(i)+'.dat'
        f = open(os.path.join(output_filename,sub_filename),'w')
        files.append(f)
    
    sub_num = len(dat)/total_subs
    count = 0
    j=0
    for i in range(0,len(dat)):
        ra = dat['RA'][i]
        dec = dat['DEC'][i]
        thi = dec*pi/180.
        phi = ra*pi/180.
        files[j].write(str(sin(thi))+' '+str(cos(thi))+' '+str(sin(phi))+' '+str(cos(phi))+' 1\n')
        #files[j].write(str(thi)+' '+str(phi)+'\n')
        if j != total_subs-1:
            count+=1
            if dat['pixl'][i] > pixflag[j] :
                j+=1
                print count
                count=0
    for i in range(0,total_subs):
        files[i].close()

survey_name = surveynames()
survey_name.eboss_step5()
pix_flag = SprtFile_nest(survey_name.uniform,20)
SprtFile_nest_dat(survey_name.data,20,pix_flag)
