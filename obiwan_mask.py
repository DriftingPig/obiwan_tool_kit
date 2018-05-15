import astropy.io.fits as fits
import os
import glob
import numpy as np
import pymangle
import subprocess
#input files
class survey():
    def __init__(self,surveyname):#eg. surveynames.eboss()
        self.obiwan = surveyname.obiwan
        self.uniform = surveyname.uniform
        self.data = surveyname.data
	self.dir = surveyname.dir
    def get_obiwan(self):
        hdu = fits.open(self.obiwan)
	dat = hdu[1].data
	hdu.close()
        ObiwanMask = dat['obiwan_mask']
        Selection = (ObiwanMask%128==1)
        obiwan_dat = np.array(dat[Selection])
        return obiwan_dat
    def get_data(self):
	hdu = fits.open(self.data)
	dat=hdu[1].data
	hdu.close()
        return dat
    def get_uniform(self):
	hdu = fits.open(self.uniform)
	dat=hdu[1].data
	hdu.close() 
        return dat

class surveynames():
    def __init__(self):
        self.dir = ''
        self.obiwan = ''
        self.uniform = ''
        self.data = ''
        self.comment = ''
	self.star_map = '/global/homes/h/huikong/eboss/LSSanalysis/maps/allstars17.519.9Healpixall256.dat'
	self.ext_map = '/global/homes/h/huikong/eboss/LSSanalysis/maps/healSFD_r_256_fullsky.dat'
	self.anand_map = '/global/homes/h/huikong/eboss/LSSanalysis/maps/ELG_hpsyst.nside256.fits'
    def eboss(self):
        self.dir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/randoms_eboss_041118/'
        self.obiwan = self.dir+'randoms_subset.fits'
        self.uniform = self.dir+'randoms_subset.fits'
        self.data = self.dir+'eBOSS_ELG_full_ALL_v1_1.dat.fits'
    def eboss_step2(self):
	self.dir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/randoms_eboss_041118/'
	self.obiwan = self.dir+'obiwan_masked_050718.fits'
	self.uniform = self.dir+'uniform_masked_050718.fits'
	self.data = self.dir+'data_masked_050718.fits'
	self.comment = 'obiwan_mask.py: vetomask and footprint mask (for eboss23)'
    def eboss_step3(self):
	self.dir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/randoms_eboss_041118/'
	self.obiwan = self.dir+'obiwan_step3.fits'
	self.uniform = self.dir+'uniform_step3.fits'
	self.data = self.dir+'data_step3.fits'
	self.comment='obiwan_bricks.py: cut the bricks to make three files consistent'
    def eboss_step4(self):
	self.dir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/randoms_eboss_041118/'
	self.obiwan = self.dir+'obiwan_step4.fits'
	self.uniform = self.dir+'uniform_step4.fits'
	self.data = self.dir+'data_step4.fits'
	self.comment='after elg selection, processed in file ELG_selection.py'
    def eboss_step5(self):
	self.dir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/randoms_eboss_041118/'
	self.obiwan = self.dir+'obiwan_step4.fits'
	self.uniform = self.dir+'uniform_step5.fits'
	self.data = self.dir+'data_step4.fits'
	self.comment='write weights to the uniform randoms'
#mask
class masks():
    def __init__(self):
        self.maskdir = ''
        self.vetomask = []
        self.footmask = []
    def eboss23(self):
        self.maskdir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/mask/eboss/'
        self.init_eboss_vetomask()
        self.init_eboss23_mask()
    def init_eboss_vetomask(self):
        eboss_mask = glob.glob(self.maskdir+'*eboss*')
        all_mask = glob.glob(self.maskdir+'*')
        veto_mask = np.array(list(set(all_mask)-set(eboss_mask)))
        self.vetomask = veto_mask
    def init_eboss23_mask(self):
        eboss_mask23 = glob.glob(self.maskdir+'*eboss23*')
        self.footmask.append(eboss_mask23[0])
        o
def mask(inl,maskl,md='veto',upper = True):
        if upper:
            RA = 'RA'
            DEC = 'DEC'
        else:
            RA = 'ra'
            DEC = 'dec'
            
        if md == 'foot':
                keep = np.zeros(inl.size,dtype='bool')  #object is outside of footprint unless it is found
        if md == 'veto':
                keep = np.ones(inl.size,dtype='bool')   #object is outside of veto mask unless it is found

        for mask in maskl:
                mng = pymangle.Mangle(mask)
                polyid = mng.polyid(inl[RA],inl[DEC])
                if md == 'foot':
                        keep[polyid!=-1] = True #keep the object if a polyid is found
                if md == 'veto':
                        keep[polyid!=-1] = False #do not keep the object if a polyid is found   
                print(mask+' done')
        return keep   
    

#mask the data    
def mask_survey(SURVEYNAME,TYPE,output_name,upper = True):
    if SURVEYNAME == 'eboss23':
        eboss = surveynames()
        eboss.eboss()
        eboss_survey = survey(eboss)
        if TYPE == 'obiwan':
            survey_dat = eboss_survey.get_obiwan()
        elif TYPE == 'uniform':
            survey_dat = eboss_survey.get_uniform()
        elif TYPE == 'data':
            survey_dat = eboss_survey.get_data()
        else:
            raise ValueError('ERROR01!')
        survey_mask = masks()
        survey_mask.eboss23()
    else:
        raise ValueError('ERROR02!')
        return False
    
    
    dat_vetol = mask(survey_dat,survey_mask.vetomask,'veto',upper)
    dat_footl = mask(survey_dat,survey_mask.footmask,'foot',upper)
    col_dat_vetol = fits.Column(name='veto_mask', format='B', array = dat_vetol)
    col_dat_footl = fits.Column(name='foot_mask', format='B', array = dat_footl)
    col_dat_orig = fits.ColDefs(np.array(survey_dat))
    col_dat_mask = col_dat_orig.add_col(col_dat_vetol).add_col(col_dat_footl)
    dat_temp = fits.BinTableHDU.from_columns(col_dat_mask).writeto('./temp.fits',overwrite = True)
    dat_masked = fits.open('./temp.fits')[1].data
    subprocess.call(["rm","temp.fits"])
    mask1 = dat_masked['veto_mask']
    mask2 = dat_masked['foot_mask']
    if upper:
	DEC = 'DEC'
    else:
	DEC = 'dec'
    mask_dec = dat_masked[DEC]
    mask_sel = (mask1==True) & (mask2==True) & (mask_dec>14.05)
    dat_final = np.array(dat_masked[mask_sel])
    print(len(dat_final),len(dat_masked),len(survey_dat))
    fits.BinTableHDU.from_columns(fits.ColDefs(np.array(dat_final))).writeto(output_name,overwrite=True)
    
def eboss23_mask_main():
    out_dir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/randoms_eboss_041118/'
    mask_survey('eboss23','obiwan',out_dir+'obiwan_masked_050718.fits',upper = False)
    mask_survey('eboss23','uniform',out_dir+'uniform_masked_050718.fits',upper = False)
    mask_survey('eboss23','data',out_dir+'data_masked_050718.fits',upper = True)
    
#eboss23_mask_main()

def file_output(data_array,file_name):
	fits.BinTableHDU.from_columns(fits.ColDefs(np.array(data_array))).writeto(file_name,overwrite=True)	

def add_column(data_set,extra_col,col_name,col_format):
	col_orig = fits.ColDefs(np.array(data_set))
	col_added = fits.Column(name=col_name, format=col_format, array = extra_col)
	col_final = col_orig.add_col(col_added)
	dat_temp = fits.BinTableHDU.from_columns(col_final).writeto('./temp.fits',overwrite = True)
	dat_final = fits.open('./temp.fits')[1].data
	subprocess.call(["rm","temp.fits"])
	return dat_final

def add_column_with_saved_file(data_set,extra_col,col_name,col_format,output_name,mode):
        col_orig = fits.ColDefs(np.array(data_set))                                       
        col_added = fits.Column(name=col_name, format=col_format, array = extra_col)      
        col_final = col_orig.add_col(col_added)
        dat_temp = fits.BinTableHDU.from_columns(col_final).writeto(output_name,overwrite = True)
        hdu = fits.open(output_name,mode = mode)
        return hdu
