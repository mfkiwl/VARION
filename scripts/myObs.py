# \author Giorgio Savastano, 2015. giorgio.savastano(at)uniroma1.it
#
# -------------------------------------------------------------------------
#
# Copyright (C) 2015-2016  (see AUTHORS file for a list of contributors)
#
# VARION is a opean source software for GNSS processing
#
# This file is part of VARION.
#
# VARION is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# VARION is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with VARION. If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------
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
