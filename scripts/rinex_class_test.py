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


rinex = RinexFile('gs190700.11o')

