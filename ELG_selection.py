from obiwan_mask import *

def select_ELG( dat ,tractor=True):
        """
        Given the path to a tractor catalog, it returns two sub catalogs with the eBOSS ELG selections applied (NGC and SGC).
        """
        if tractor==True:
		prefix='tractor_'
	else:
		prefix=''
        # the color color selection
	if tractor==True:
        	g     = 22.5 - 2.5 * np.log10(dat[prefix+'flux_g'] / dat[prefix+'mw_transmission_g'])
        	r_mag = 22.5 - 2.5 * np.log10(dat[prefix+'flux_r'] / dat[prefix+'mw_transmission_r'])
        	z_mag = 22.5 - 2.5 * np.log10(dat[prefix+'flux_z'] / dat[prefix+'mw_transmission_z'])
        	# opens the tractor file
        	gr = g - r_mag
        	rz = r_mag - z_mag
	else:
		g=dat['g']
		gr=dat['gr']
		rz=dat['rz']
        color_sgc = (g>21.825)&(g<22.825)&(-0.068*rz+0.457<gr)&(gr< 0.112*rz+0.773) &(0.218*gr+0.571<rz)&(rz<-0.555*gr+1.901)
        color_ngc = (g>21.825)&(g<22.9)  &(-0.068*rz+0.457<gr)&(gr< 0.112*rz+0.773) &(0.637*gr+0.399<rz)&(rz<-0.555*gr+1.901)
        # the junk rejection criterion
	if tractor==True:
        	noJunk = (dat[prefix+'anymask_g']==0) & (dat[prefix+'anymask_r']==0) & (dat[prefix+'anymask_z']==0) #& (dat['TYCHO2INBLOB']==False)
	else:
		noJunk=True
        # the low depth region rejection
        value_g=dat[prefix+'psfdepth_g']
        value_r=dat[prefix+'psfdepth_r']
        value_z=dat[prefix+'psfdepth_z']
        gL = 62.79716079
        rL = 30.05661087
        zL_ngc = 11.0
        zL_sgc = 12.75
        depth_selection_ngc = (value_g > gL) & (value_r > rL) & (value_z > zL_ngc)
        depth_selection_sgc = (value_g > gL) & (value_r > rL) & (value_z > zL_sgc)
        # final selection boolean array :
        selection_sgc =(noJunk)&(color_sgc)&(depth_selection_sgc)
        selection_ngc =(noJunk)&(color_ngc)&(depth_selection_ngc)
        # returns the catalogs of ELGs
        if len(selection_sgc.nonzero()[0])>0 or  len(selection_ngc.nonzero()[0])>0 :
                flag = True
                return flag, dat[selection_ngc], dat[selection_sgc]
        else :
                flag = False
                return flag, dat[selection_ngc], dat[selection_sgc]

def main():
	eboss_name=surveynames()
	eboss_name.eboss_step3()
	eboss_survey=survey(eboss_name)
	obiwan=eboss_survey.get_obiwan()
	uniform=eboss_survey.get_uniform()
	data=eboss_survey.get_data()
	_,obiwan_elg,_ = select_ELG( obiwan )
	_,uniform_elg,_ = select_ELG( uniform )
	_,data_elg,_ = select_ELG( data ,tractor=False)
	file_output(obiwan_elg,eboss_survey.dir+'obiwan_step4.fits')
	file_output(uniform_elg,eboss_survey.dir+'uniform_step4.fits')
	file_output(data_elg,eboss_survey.dir+'data_step4.fits')
	return True
main()
