import numpy as np

class RinexFile:
	"""docstring for RinexFile"""
	def __init__(self, rinex):
		self.nam = rinex
		self.prg = ""
		self.int = ""
		self.typ = ""
		self.xyz = ""
	def PROGRAM_GENERATOR(self):
		'''
		Method that exstracts the name of the Program Generator --> e.g. Spider, teqc...
		'''
		with open(self.nam,'r') as f:
			for i in xrange(0,3):
				line = f.readline()
				if 'PGM' in line:
					ln = (line).split()
					prog = ln[0]			
					self.prg = prog
					break
	def COORD_XYZ(self):
		'''
		Method that exstracts the x, y, z, (WGS84) --> APPROX of the receiver
		'''
		with open(self.nam, "r") as f:
			for i in xrange(0,35):
					line = f.readline() 
					if 'APPROX' in line:
						x = float(line[1:14])
						y = float(line[15:29])
						z = float(line[29:43])
						self.xyz = np.asarray([x, y, z])
						break
	def TYPE_OBS(self):
		'''
		output: name of the Program Generator --> e.g. Spider, teqc...
		'''
		import re
		with open(self.nam, "r") as f:
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
					self.typ = np.asarray([c1_index, l1_index, l2_index])
					break

	def INTERVAL(self):   ### DEBUGGG AND TEST
		'''
		Function that returns the inteval of the obs.
		'''
		with open(self.nam, "r") as f: 
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
					self.int = interval
					break


rinex = RinexFile('gs190700_test.11o')

