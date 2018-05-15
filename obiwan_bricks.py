from obiwan_mask import *

class bricks():
	def __init__(self,survey):
		self.obiwan = survey.get_obiwan()
		self.uniform = survey.get_uniform()
		self.data = survey.get_data()
		print len(self.obiwan),len(self.uniform),len(self.data)
		self.dir = survey.dir
	def bricklist_set(self,data_type):
		print len(data_type)
		brickname = data_type['brickname']
		return set(brickname)
	def brick_overlap(self):
		obiwan_set = self.bricklist_set(self.obiwan)
		uniform_set = self.bricklist_set(self.uniform)
		data_set = self.bricklist_set(self.data)
		return obiwan_set.intersection(uniform_set).intersection(data_set)
	def brick_mask_list(self,data_type,overlap_set):
		brick_mask = np.ones(len(data_type))
		for i in range(len(data_type)):
			if overlap_set-set([data_type['brickname'][i]]) == overlap_set:
				brick_mask[i]=0
		return brick_mask
	def add_brick_mask_col(self,data_type,brick_mask_list):
		from astropy.table import Column
		from astropy.table import Table
		col1 = fits.ColDefs(np.array(data_type))
		col2 = fits.Column(name='brick_mask', format='B', array = brick_mask_list)
		cols = col1.add_col(col2)
		fits.BinTableHDU.from_columns(cols).writeto('./temp.fits',overwrite = True)
		dat_masked = fits.open('./temp.fits')[1].data
		subprocess.call(["rm","temp.fits"])
		mask = dat_masked['brick_mask']
		mask_sel = (mask==True)
		dat_final = np.array(dat_masked[mask_sel])
		return dat_final
	def eboss_brick_analysis(self):
		overlapped_brick = self.brick_overlap()
		print len(overlapped_brick)
		obiwan_brick_mask_list = self.brick_mask_list(self.obiwan,overlapped_brick)
		obiwan_brick_mask_col = self.add_brick_mask_col(self.obiwan,obiwan_brick_mask_list)
		fits.BinTableHDU.from_columns(fits.ColDefs(np.array(obiwan_brick_mask_col))).writeto(self.dir+'obiwan_step3.fits',overwrite = True)

		uniform_brick_mask_list = self.brick_mask_list(self.uniform,overlapped_brick)
		uniform_brick_mask_col = self.add_brick_mask_col(self.uniform,uniform_brick_mask_list)
		fits.BinTableHDU.from_columns(fits.ColDefs(np.array(uniform_brick_mask_col))).writeto(self.dir+'uniform_step3.fits',overwrite = True)

		data_brick_mask_list = self.brick_mask_list(self.data,overlapped_brick)
		data_brick_mask_col = self.add_brick_mask_col(self.data,data_brick_mask_list)
		fits.BinTableHDU.from_columns(fits.ColDefs(np.array(data_brick_mask_col))).writeto(self.dir+'data_step3.fits',overwrite = True)

def main():
	#get files
	eboss_step2 = surveynames()
	eboss_step2.eboss_step2()
	eboss_survey=survey(eboss_step2)
	eboss_bricks=bricks(eboss_survey)
	eboss_bricks.eboss_brick_analysis()

main()
	
