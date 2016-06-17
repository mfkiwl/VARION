import numpy as np

mask_var = data_varion['prn'] == '10'
mask_vad = data['prn'] == 'G10'
data_varion = np.genfromtxt( 'ainp3020.12n_VARION.sky', skip_header=1, usecols=(0,1,2,3,4,5), dtype=[('prn','S5'),('sow','f8'),('x','f8'),('y','f8'),('z','f8'),('toe','i8')] )
data = np.genfromtxt( 'ainp3020.12o_L3_G_L1W+L2W_ALL.Orb_G', skip_header=1, usecols=(0,1,2,3,4,6), dtype=[('prn','S5'),('sow','f8'),('x','f8'),('y','f8'),('z','f8'),('toe','i8')] )
