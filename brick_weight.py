from obiwan_mask import *
def main():
	#read files
	eboss_name=surveynames()
	eboss_name.eboss_step4()
	eboss_survey=survey(eboss_name)
	obiwan = eboss_survey.get_obiwan()
	uniform = eboss_survey.get_uniform()
	data = eboss_survey.get_data()
	#brick process
	bricklist = list(set(obiwan['brickname']))#list of bricks
	brick_length = len(bricklist)#length of bricks
	#obiwan_average = float(len(obiwan))/float(brick_length)
	#uniform_average = float(len(uniform))/float(brick_length)
	#brick loop, weight per brick
	obiwan_brickname = obiwan['brickname']
	uniform_brickname = uniform['brickname']
	weight_list = []
	for i in range(brick_length):
		SEL_obiwan = (obiwan_brickname==bricklist[i])
		SEL_uniform = (uniform_brickname==bricklist[i])
		brick_array_obiwan = np.array(obiwan[SEL_obiwan])
		brick_array_uniform = np.array(uniform[SEL_uniform])
		per_brick_length_obiwan = len(brick_array_obiwan)
		per_brick_length_uniform = len(brick_array_uniform)
		brick_weight = float(per_brick_length_obiwan)/float(per_brick_length_uniform)
		weight_list.append(brick_weight)
	#test
	'''
	from matplotlib import pyplot as plt
	plt.hist(weight_list,bins=100)
	plt.show()
	'''
	weight_init = np.ones(len(uniform))
	add1 = add_column(uniform,weight_init,'obiwan_weight', 'D')
	ids = np.arange(len(uniform))
	hdu_uniform_col_added = add_column_with_saved_file(add1,ids,'ids', 'I',eboss_survey.dir+'uniform_step5.fits','update')
	dat_uniform_col_added = hdu_uniform_col_added[1].data
	w_brickname = dat_uniform_col_added['brickname']
	
	for i in range(brick_length):
		SEL_brick = (w_brickname==bricklist[i])
		selected = np.array(dat_uniform_col_added[SEL_brick])
		for j in range(len(selected)):
			index = selected['ids'][j]
			dat_uniform_col_added['obiwan_weight'][index]=weight_list[i]
	hdu_uniform_col_added.flush()
	#test
	'''
	from matplotlib import pyplot as plt
	plt.hist(dat_uniform_col_added['obiwan_weight'])
	plt.show()
	'''
main()
