##########

class Observation(object):
    """A class that makes various tasty fruits."""
    def __init__(self, sat, epoch, freq, value):
        self.sat    = sat
        self.epoch  = epoch
        self.freq   = freq
        self.value  = value

    def description(self):
        print "Satellite: %s at epoch: %s, signal: %s." % (self.sat, self.epoch, self.freq)
        
##########

def geometry_free(L1,L2):
	a =  1
	b = -1
	L4 = a*L1 + b*L2
	return L4

##########

def obs_sat(sats_arr, ora_arr, sod_arr, l1_arr, l2_arr, name_sat):
    """
    input:
        - array of all satellites for all epochs
        - array of all L1 obs
        - array of all L1 obs
        - mask_sat = witch satellites have to be selected
    output:
        - array of the VAD geometry-free for the selected sat
        - array of the hour of each epoch
        - arrayof the second of the day
    """
    mask_sat = (sats_arr == name_sat) 
    sat    = sats_arr[mask_sat]
    ora    = ora_arr[mask_sat]
    sod    = sod_arr[mask_sat]
    sat_l1 = l1_arr[mask_sat]
    sat_l2 = l2_arr[mask_sat]
    sat_l4 = geometry_free(sat_l1,sat_l2)
    sat_VAD= sat_l4[1:len(sat_l4)] - sat_l4[0:len(sat_l4)-1]  
    return sat_VAD, ora[1:len(ora)], sod[1:len(sod)]
