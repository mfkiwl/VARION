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
import numpy as np
##
class RinexFile:
	"""
	Class for RINEX attributes. 
	The class contins the following methods automatically called by the constructor:
			- PROGRAM_GENERATOR()  -->  set the self.prg with the Program Generator name 
			- COORD_XYX()          -->  set the self.xyz with the APPROX coord of the receiver
			- TYPE_OBS()           -->  set the self.typ with the order of the obs (C1, L1, L2 ...)
			- INTERVAL()           -->  set the self.int with the time interval of the obs
	the class also contains the following method:
			- READ_RINEX()         -->  set the self.data with al the obs in the RINEX file
	"""
	# Init
	def __init__(self, rinex):
		self.nam = rinex
		# Init properties of the object by calling some methods
		def PROGRAM_GENERATOR(file_nam):
			'''
			Method that exstracts the name of the Program Generator --> e.g. Spider, teqc...
			'''
			with open(file_nam,'r') as f:
				for i in xrange(0,3):
					line = f.readline()
					if 'PGM' in line:
						ln = (line).split()
						prog = ln[0]			
						break
			return prog
		def COORD_XYZ(file_nam):
			'''
			Method that returns the x, y, z, (WGS84) --> APPROX of the receiver
			'''
			with open(file_nam, "r") as f:
				for i in xrange(0,35):
						line = f.readline() 
						if 'APPROX' in line:
							x = float(line[1:14])
							y = float(line[15:29])
							z = float(line[29:43])
							break
			return np.asarray([x, y, z])
		def TYPE_OBS(file_nam):
			'''
			Method that returns the name of the Program Generator --> e.g. Spider, teqc...
			'''
			import re
			with open(file_nam, "r") as f:
				for i in xrange(35):
					lns = f.readline()
					if "TYPES OF OBSERV" in lns:
						obs_ord = re.findall('\S+\S+',lns)
						obs_orb_arr = np.array(obs_ord)
						c1_mas = (obs_orb_arr == "C1")
						l1_mas = (obs_orb_arr == "L1")
						l2_mas = (obs_orb_arr == "L2")
						index = np.arange(len(obs_orb_arr))
						c1_index = index[c1_mas][0]
						l1_index = index[l1_mas][0]
						l2_index = index[l2_mas][0]
						break
			return np.asarray([c1_index, l1_index, l2_index]), lns[4:6]
		### DEBUGG and TEST
		def GPSTIME(file_nam):
			'''
			'''
			import re
			import GPS
			with open(file_nam, "r") as f: 
				for i in range(100):
					lns = f.readline()
					if "TIME OF FIRST OBS" in lns:
						lns = lns.split()
						year  = int(lns[0])
						month = int(lns[1])
						day   = int(lns[2])
						gps_week_ref, secsOfWeek, gpsDay, gpsSOD = GPS.gpsFromUTC(year, month, day,  0,  0,  0, leapSecs=0)
						break
			return gps_week_ref, secsOfWeek
		def INTERVAL(file_nam):   ### DEBUGGG AND TEST
			'''
			Function that returns the inteval of the obs
			'''
			import re
			with open(file_nam, "r") as f: 
				while True:
					lns = f.readline()
					if "INTERVAL" in lns:
						inter = re.findall('\S+\.\S+',lns)
						interval = float(inter[0]) 
						self.int = interval
						break
					elif "END OF HEADER" in lns:
						sod = []
						for i in xrange(100):
							lns = f.readline()
							if "COMMENT" in lns:
								continue 
							else:
								pass  
							# START stocking the obs                           
							ln = lns.split()                         # split the line[i]
							if len(ln) >= 8:                         # check the number of the spli
								sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))  
						sod_array = np.asarray(sod)             
						interval_array = sod_array[1:len(sod_array)] - sod_array[0:len(sod_array)-1]   
						interval = np.median( interval_array )         # use the median to define the interval in order to remove the outlayers
						break
			return interval
		#
		self.prg = PROGRAM_GENERATOR( self.nam )
		self.xyz = COORD_XYZ( self.nam )
		self.int = INTERVAL( self.nam )
		self.typ, self.obs = TYPE_OBS( self.nam )
		self.gps_week_ref, self.gps_sow_ref = GPSTIME( self.nam )
		self.data = ""
	#######################################################################
	def READ_RINEX(self):
		'''
		Function to stock the obs of the rinex
		'''
		import re
		import datetime
		L1 = 1.57542e9    # Hz
		L2 = 1.22760e9    # Hz
		c  = 299792458    # m/s
		#
		lam1=c/L1         # m
		lam2=c/L2         # m
		##
		num_lin = 1
		if int(self.obs)<=5:
				num_lin = 2	
		#
		f   = open(self.nam, "r")
		lns = f.readlines()
		f.close()
		#
		sats = []; ora  = []; sod  = []; c1   = []; l1   = []; l2   = []
		s_sats= 0; count = 0
		# Skip header
		while "END OF HEADER" not in lns[count]:
				count += 1
		count += 1
		for i in xrange(count,len(lns)): 
					## skip the comment inside the file
					if "COMMENT" in lns[i]: continue 
					# START stocking the obs                           
					ln = lns[i].split()
					## DEBUG  ---> per quei file osservazioni che hanno buchi sui satelliti con lo zero --- >= 8
					if len(ln) >= 8:
							new_ln = ln[7] 
							for yy in xrange( 1, (len(ln)-8)+1):             # se len(ln) e' proprio uguale a 8 non dovrebbe fare nulla
								if self.prg == 'DAT2RIN':
									if len(ln[7+yy]) == 1:
										new_ln += "G0" + ln[7+yy]
									else:
										new_ln += "G" + ln[7+yy]
								else:
										new_ln += "0" + ln[7+yy] 
							ln[-1] = new_ln    
						## count how many S sats there are
							if "S" in ln[-1]:
								s_sats=0
								for s in ln[-1]:
									if s == "S":
										s_sats+=1
									else: pass
							##
							if (len(ln[-1]) - 1) % 3 == 0:                                         
									n_sats = int(ln[-1][0])                                         
									skip_num = 1
							else:                             
									n_sats = int(ln[-1][0:2])
									# skip the first two characters because the sat number is > 10 
									skip_num = 2
							##
							if n_sats <= 12:
									# num_lin to read also RINEX SPIDER
									lin_epoc = (n_sats-s_sats)*(2/num_lin)    
									for k in xrange(1, lin_epoc+1):    
											# salto le righe pari
											if num_lin != 2:
												if k % 2 == 0: continue   
											# use regular expression to select only the float number
											obs = re.findall('\S+\.\S+', lns[i+k])
											if len(obs) <= 3: continue
											# split the obs and truncate the number to the 3rd decimal (no rounding)
											c1_before_dec, c1_after_dec = str(obs[ self.typ[0] ]).split('.')                             
											l1_before_dec, l1_after_dec = str(obs[ self.typ[1] ]).split('.')     
											l2_before_dec, l2_after_dec = str(obs[ self.typ[2] ]).split('.')
											c1.append(   float('.'.join(( c1_before_dec, c1_after_dec[0:3])))   )    
											l1.append(   float('.'.join(( l1_before_dec, l1_after_dec[0:3])))  *lam1   )
											l2.append(   float('.'.join(( l2_before_dec, l2_after_dec[0:3])))  *lam2   )            
											# slicing
											if num_lin == 1:
											    sats.append( ln[-1][skip_num:len(ln[-1])][(k)/(2/num_lin)*3:(k)/(2/num_lin)*3+3] )      
											else:
											    sats.append( ln[-1][skip_num:len(ln[-1])][(k-1)*3:(k-1)*3+3] )
											timing = (  datetime.time( int(ln[3]), int(ln[4]) , int(float(ln[5][0:2])) ) )
											ora.append ( timing.isoformat() )
											sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))
							else:                                                      # Check the number of sats, if >12 use continuation lines
									ln_2 = lns[i+1].split()                     # splitto anche la riga dopo
									lin_epoc = (n_sats-s_sats)*(2/num_lin)
									for k in xrange(1, lin_epoc+1):                      # if n_sats=12, than k = [1,3,5,...,23]
											if num_lin != 2:
												if k % 2 == 0: continue 
											obs = re.findall('\S+\.\S+',lns[i+1+k])
											if len( obs ) <= 3: continue
											c1_before_dec, c1_after_dec = str(obs[ self.typ[0] ]).split('.')                         
											l1_before_dec, l1_after_dec = str(obs[ self.typ[0] ]).split('.')     
											l2_before_dec, l2_after_dec = str(obs[ self.typ[0] ]).split('.')
											c1.append(   float('.'.join(( c1_before_dec, c1_after_dec[0:3])))   )    
											l1.append(   float('.'.join(( l1_before_dec, l1_after_dec[0:3])))  *lam1   )
											l2.append(   float('.'.join(( l2_before_dec, l2_after_dec[0:3])))  *lam2   )                     
											timing = (  datetime.time( int(ln[3]), int(ln[4]) , int(float(ln[5][0:2])) ) )
											ora.append ( timing.isoformat() )
											sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))     
											if k <= 23:
											            sats.append(ln[-1][2:len(ln[-1])][(k)/(2/num_lin)*3:(k)/(2/num_lin)*3+3])  
											else:
											            sats.append(ln_2[0][((k+1)/(2/num_lin)-13)*3:((k+1)/(2/num_lin)-13)*3+3])
		sats_arr = np.asarray(sats)                                          
		ora_arr  = np.asarray(ora)
		sod_arr  = np.asarray(sod)
		c1_arr = np.asarray(c1)
		l1_arr = np.asarray(l1)
		l2_arr = np.asarray(l2)
		self.data = sats_arr, ora_arr, sod_arr, l1_arr, l2_arr, c1_arr,
