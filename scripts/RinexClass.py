# -*- coding: utf-8 -*-
'''


created  by  Giorgio Savastano in 2015    <giorgio.savastano(at)uniroma1.it>
modified by  Michela Ravanelli in 2018    <michela.ravanelli(at)uniroma1.it>

 -------------------------------------------------------------------------

 Copyright (C) 2015-2020  (see AUTHORS file for a list of contributors)

 VARION is a open source software for GNSS processing

 This file is part of VARION.

 VARION is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 VARION is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with VARION. If not, see <http://www.gnu.org/licenses/>.

 -------------------------------------------------------------------------
'''
 
import numpy as np
import pdb
import pandas as pd
import time
import collections
from decimal import Decimal


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
		def RINEX_VERSION(file_nam):
			'''
			Method that exstracts the name of the Program Generator --> e.g. Spider, teqc...
			'''
			with open(file_nam,'r') as f:
				line = f.readline()
				#pdb.set_trace()
				if 'RINEX VERSION' in line:
						ln = (line).split()
						version = ln[0]			

			return version
		def PROGRAM_GENERATOR(file_nam):
			'''
			Method that exstracts the name of the Program Generator --> e.g. Spider, teqc...
			'''
			with open(file_nam,'r') as f:
				for i in xrange(0,6):
					line = f.readline()
					if 'PGM' in line:
						ln = (line).split()
						prog = ln[0]			
						break
					else:
						prog='no prog generator'
			return prog
		def COORD_XYZ(file_nam):
			'''
			Method that returns the x, y, z, (WGS84) --> APPROX of the receiver
			'''
			with open(file_nam, "r") as f:
				for i in xrange(0,40):
						line = f.readline() 
						if 'APPROX' in line:
							x = float(line[1:14])
							y = float(line[15:29])
							z = float(line[29:43])
							break
						else:
							x = float('Nan')
							y = float('Nan')
							z = float('Nan')
			return np.asarray([x, y, z])
		def TYPE_OBS( file ):
			'''
			Method that returns the name of the Program Generator --> e.g. Spider, teqc...
			'''
			import re
			with open( file.nam, "r") as f:
				lines=f.readlines()

				obs_ord_gps = np.zeros(10) #versione 3 questo serve per tutelarmi laddove non ho qualche osservazione 
				obs_ord_sbas = np.zeros(10)
				obs_ord_china = np.zeros(10)
				obs_ord_gali=np.zeros(10)

				obs_v2_number=np.zeros(2) #versione2
				obs_ord_v2=np.zeros(10)


				for i in xrange(260):
					lns = f.readline()
					#print i, lns 

					if float(file.vrs) < 3:

						if len(lines[i].split())<3:
							continue
						else: 
							if lines[i].split()[-3]=='TYPES' and lines[i].split()[-2]=='OF' and lines[i].split()[-1]=='OBSERV' :
								obs_v2_number=lines[i].split()[0]

								obs_ord = lines[i].split()
								obs_orb_arr = np.array(obs_ord)
								try:


										if Decimal(float(obs_v2_number)/9)%1!=0:
											righe=int(obs_v2_number)/9								
										else:
											righe=int(obs_v2_number)/9-1


										obs_ord1=[]

										for k in xrange(0,righe+1):
											obs_ord1.append(lines[i+k].split())
										obs_orb_arr=np.asarray(np.hstack(obs_ord1))

										c1_mas = (obs_orb_arr == "C1") 
										if len(obs_orb_arr[c1_mas])==0:
											c1_mas = (obs_orb_arr == "P1")
										l1_mas = (obs_orb_arr == "L1")
										

										a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='OF')[0])
										b=np.delete(a,np.where(a=='/')[0])
										c=np.delete(b,np.where(b=='OBSERV')[0])
										d=np.delete(c,np.where(c=='TYPES')[0])
										obs_ord_v2=np.delete(d,np.where(d=='#')[0])
	

										
								except ValueError:
									continue
						if lines[i].split()[0]=='END' and lines[i].split()[1]=='OF' and lines[i].split()[2]=='HEADER' :
							break


					else:



						obs_ord = lines[i].split()
						obs_orb_arr = np.array(obs_ord)
						

						#print i, lines[i], obs_ord_gps, obs_ord_sbas,obs_ord_china

						if lines[i].split()[0]=='G' and lines[i].split()[-1]!='SHIFT':
							obs_number=lines[i].split()[1]


							if int(obs_number)/13!=0:
								righe=int(obs_number)/13
								obs_ord1=[]
								for k in xrange(0,righe+1):
									obs_ord1.append(lines[i+k].split())
								
								obs_orb_arr=np.asarray(np.hstack(obs_ord1))
								c1_mas = (obs_orb_arr == "C1C")
								l1_mas = (obs_orb_arr == "L1C")
								l2_mas = (obs_orb_arr == "L2W")
								c2_mas = (obs_orb_arr == "C2W")
								index = np.arange(len(obs_orb_arr))
								c1_index = index[c1_mas][0] -1
								l1_index = index[l1_mas][0] -1
								l2_index = index[l2_mas][0] -1
								c2_index = index[c2_mas][0] -1

								a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
								b=np.delete(a,np.where(a=='/')[0])
								c=np.delete(b,np.where(b=='OBS')[0])
								d=np.delete(c,np.where(c=='TYPES')[0])
								obs_ord_gps=np.delete(d,np.where(d=='#')[0])




							else:

								c1_mas = (obs_orb_arr == "C1C") 
								l1_mas = (obs_orb_arr == "L1C") 
								l2_mas = (obs_orb_arr == "L2W") 
								c2_mas = (obs_orb_arr == "C2W")
								index = np.arange(len(obs_orb_arr))
								c1_index = index[c1_mas][0] -1
								l1_index = index[l1_mas][0] -1
								try:
									l2_index = index[l2_mas][0] -1
									c2_index = index[c2_mas][0] -1
								except IndexError:
									l2_index = float('Nan')
									c2_index = float('Nan')
								
								a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
								b=np.delete(a,np.where(a=='/')[0])
								c=np.delete(b,np.where(b=='OBS')[0])
								d=np.delete(c,np.where(c=='TYPES')[0])
								obs_ord_gps=np.delete(d,np.where(d=='#')[0])

								

						elif lines[i].split()[0]=='S' and  lines[i].split()[-1]!='SHIFT' :
								obs_number=lines[i].split()[1]
								if int(obs_number)/13!=0:
									righe=int(obs_number)/13
									obs_ord1=[]
									for k in xrange(0,righe+1):
										obs_ord1.append(lines[i+k].split())
									obs_orb_arr=np.asarray(np.hstack(obs_ord1))
									c1_mas = (obs_orb_arr == "C1C")
									l1_mas = (obs_orb_arr == "L1C")
									c5_mas = (obs_orb_arr == "C5I")
									l5_mas = (obs_orb_arr == "L5I")
									index = np.arange(len(obs_orb_arr))
									c1_S_index = index[c1_mas][0]
									l1_S_index = index[l1_mas][0]
									try:
											c5_S_index = index[c5_mas][0] -1
											l5_S_index = index[l5_mas][0] -1
									except IndexError:
										c5_S_index = float('Nan')
										l5_S_index = float('Nan')

									a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
									b=np.delete(a,np.where(a=='/')[0])
									c=np.delete(b,np.where(b=='OBS')[0])
									d=np.delete(c,np.where(c=='TYPES')[0])
									obs_ord_sbas=np.delete(d,np.where(d=='#')[0])

								else:
									c1_mas = (obs_orb_arr == "C1C")
									l1_mas = (obs_orb_arr == "L1C")
									c5_mas = (obs_orb_arr == "C5I")
									l5_mas = (obs_orb_arr == "L5I")
									index = np.arange(len(obs_orb_arr))
									c1_S_index = index[c1_mas][0]
									l1_S_index = index[l1_mas][0]
									try:
											c5_S_index = index[c5_mas][0] -1
											l5_S_index = index[l5_mas][0] -1
									except IndexError:
										c5_S_index = float('Nan')
										l5_S_index = float('Nan')

									a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
									b=np.delete(a,np.where(a=='/')[0])
									c=np.delete(b,np.where(b=='OBS')[0])
									d=np.delete(c,np.where(c=='TYPES')[0])
									obs_ord_sbas=np.delete(d,np.where(d=='#')[0])



						elif lines[i].split()[0]=='E' and lines[i].split()[-1]!='SHIFT':
							obs_number=lines[i].split()[1]
							if int(obs_number)/13!=0:
								righe=int(obs_number)/13
								obs_ord1=[]
								for k in xrange(0,righe+1):
									obs_ord1.append(lines[i+k].split())
								
								obs_orb_arr=np.asarray(np.hstack(obs_ord1))							
								c1_mas = (obs_orb_arr == "C1C") #puoi cambiare o aggiungere C1I se serve
								l1_mas = (obs_orb_arr == "L1C")
								l5_mas = (obs_orb_arr == "L5X")
								c5_mas = (obs_orb_arr == "C5X")
								index = np.arange(len(obs_orb_arr))
								
								try:
										c1_index = index[c1_mas][0] -1
										l1_index = index[l1_mas][0] -1
								except IndexError:

									c1_index = float('Nan')
									l1_index = float('Nan')   
								
								try:
										l5_index = index[l5_mas][0] -1
										c5_index = index[c5_mas][0] -1
								except IndexError:
									c5_index = float('Nan')
									l5_index = float('Nan')


								a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
								b=np.delete(a,np.where(a=='/')[0])
								c=np.delete(b,np.where(b=='OBS')[0])
								d=np.delete(c,np.where(c=='TYPES')[0])
								obs_ord_gali=np.delete(d,np.where(d=='#')[0])
							else:
								
								c1_mas = (obs_orb_arr == "C1C")
								l1_mas = (obs_orb_arr == "L1C")
								l5_mas = (obs_orb_arr == "L5X")
								c5_mas = (obs_orb_arr == "C5X")
								index = np.arange(len(obs_orb_arr))
								
								try:
										c1_index = index[c1_mas][0] -1
										l1_index = index[l1_mas][0] -1
								except IndexError:
										c1_index = float('Nan')
										l1_index = float('Nan')        
									
								try:		
										l5_index = index[l5_mas][0] -1
										c5_index = index[c5_mas][0] -1
									
								except IndexError:		

									c5_index = float('Nan')
									l5_index = float('Nan')

								a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
								b=np.delete(a,np.where(a=='/')[0])
								c=np.delete(b,np.where(b=='OBS')[0])
								d=np.delete(c,np.where(c=='TYPES')[0])
								obs_ord_gali=np.delete(d,np.where(d=='#')[0])




						elif lines[i].split()[0]=='C'and lines[i].split()[-1]!='SHIFT':
							obs_number=lines[i].split()[1]
							if int(obs_number)/13!=0:
								righe=int(obs_number)/13
								obs_ord1=[]
								for k in xrange(0,righe+1):
									obs_ord1.append(lines[i+k].split())
								
								obs_orb_arr=np.asarray(np.hstack(obs_ord1))
								c1_mas = (obs_orb_arr == "C1I")
								l1_mas = (obs_orb_arr == "L1I")
								if len( obs_orb_arr[obs_orb_arr == "C1I"])!=0: #nei rinex 3.01 la prima freq e' C2I, nei 3.02 e'C1I
									c1_mas = (obs_orb_arr == "C1I")
									l1_mas = (obs_orb_arr == "L1I")
								elif len( obs_orb_arr[obs_orb_arr == "C2I"])!=0:
									c1_mas = (obs_orb_arr == "C2I")
									l1_mas = (obs_orb_arr == "L2I")
								elif len( obs_orb_arr[obs_orb_arr == "C1I"])==0:
									if len( obs_orb_arr[obs_orb_arr == "C1X"])==0:
										c1_mas = (obs_orb_arr == "C1C")
										l1_mas = (obs_orb_arr == "L1C")
									else:
										c1_mas = (obs_orb_arr == "C1X")
										l1_mas = (obs_orb_arr == "L1X")
								l7_mas = (obs_orb_arr == "L7I") 
								c7_mas = (obs_orb_arr == "C7I")
								l6_mas = (obs_orb_arr == "L6I") 
								c6_mas = (obs_orb_arr == "C6I")

								index = np.arange(len(obs_orb_arr))
								c1_C_index = index[c1_mas][0] -1
								l1_C_index = index[l1_mas][0] -1

								try:
									l7_C_index = index[l7_mas][0] -1
									c7_C_index = index[c7_mas][0] -1
									l6_C_index = index[l6_mas][0] -1
									c6_C_index = index[c6_mas][0] -1
								except IndexError:
									l7_C_index = float('Nan')
									c7_C_index = float('Nan')
									l6_C_index = float('Nan')
									c6_C_index = float('Nan')

								a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
								b=np.delete(a,np.where(a=='/')[0])
								c=np.delete(b,np.where(b=='OBS')[0])
								d=np.delete(c,np.where(c=='TYPES')[0])
								obs_ord_china=np.delete(d,np.where(d=='#')[0])


							else:

								if len( obs_orb_arr[obs_orb_arr == "C1I"])!=0:
										c1_mas = (obs_orb_arr == "C1I")
										l1_mas = (obs_orb_arr == "L1I")

								elif len( obs_orb_arr[obs_orb_arr == "C2I"])!=0:
										c1_mas = (obs_orb_arr == "C2I")
										l1_mas = (obs_orb_arr == "L2I")

								elif len( obs_orb_arr[obs_orb_arr == "C1I"])==0:
									if len( obs_orb_arr[obs_orb_arr == "C1X"])==0:
										c1_mas = (obs_orb_arr == "C1C")
										l1_mas = (obs_orb_arr == "L1C")
									else:
										c1_mas = (obs_orb_arr == "C1X")
										l1_mas = (obs_orb_arr == "L1X")

								
								l7_mas = (obs_orb_arr == "L7I") 
								c7_mas = (obs_orb_arr == "C7I")
								l6_mas = (obs_orb_arr == "L6I") 
								c6_mas = (obs_orb_arr == "C6I")

								index = np.arange(len(obs_orb_arr))
								c1_C_index = index[c1_mas][0] -1
								l1_C_index = index[l1_mas][0] -1

								try:
											l7_C_index = index[l7_mas][0] -1
											c7_C_index = index[c7_mas][0] -1
											l6_C_index = index[l6_mas][0] -1
											c6_C_index = index[c6_mas][0] -1

								except IndexError:
										l7_C_index = float('Nan')
										c7_C_index = float('Nan')
										l6_C_index = float('Nan')
										c6_C_index = float('Nan')

								a=np.delete(obs_orb_arr,np.where(obs_orb_arr=='SYS')[0])
								b=np.delete(a,np.where(a=='/')[0])
								c=np.delete(b,np.where(b=='OBS')[0])
								d=np.delete(c,np.where(c=='TYPES')[0])
								obs_ord_china=np.delete(d,np.where(d=='#')[0])



					if lines[i].split()[0]=='END' and lines[i].split()[1]=='OF' and lines[i].split()[2]=='HEADER' :
						break
					#elif lines[i].split()[-1]=='SHIFT':
					#	break
					else:
						pass


			return obs_ord_gps, obs_ord_sbas,obs_ord_china, obs_ord_gali ,obs_ord_v2 
			#obs_ord_gps, obs_ord_sbas,obs_ord_china, obs_ord_gali  mi restituisce le osservazioni presenti nel rinex  per gps, sbas, beidou e galileo


		### DEBUGG and TEST
		def GPSTIME(file_nam):
			'''
			'''
			import re
			import GPS
			with open(file_nam, "r") as f: 
				for i in range(240):
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
						interval = float(lns.split()[0])
						self.int = interval						
						break
					elif "END OF HEADER" in lns:
						sod = []
						for i in xrange(120):
							lns = f.readline()

							if "COMMENT" in lns:
								continue 
							else:
								pass  
							# START stocking the obs                           
							ln = lns.split() # split the line[i]
							#print i,ln
							if float(self.vrs)<3:#rinex 2
								if len(ln)>=8:
								        if len(ln[0])==2:
									   sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3])) #rinex2
							else: #rinex3
								if len(ln)==9:
									if ln[0]== 'R' or ln[0]=='G' or ln[0]=='S' or ln[0]=='J' or ln[0]=='C' or ln[0][0]== 'R' or ln[0][0]== 'G' or ln[0]=='E' or ln[0][0]=='E' or ln[0][0]=='J'or ln[0][0]=='C' or ln[0][0]=='S' :
										continue
									sod.append(int(ln[4])*60*60 + int(ln[5])*60 + float(ln[6][0:3])) #rinex3
							if len(sod)==2:
									break
						sod_array = np.asarray(sod) 
						interval_array = sod_array[1:len(sod_array)] - sod_array[0:len(sod_array)-1]   
						interval = np.median( interval_array )         # use the median to define the interval in order to remove the outlayers
						break
			return interval



		self.vrs = RINEX_VERSION( self.nam )
		self.prg = PROGRAM_GENERATOR( self.nam )
		self.xyz = COORD_XYZ( self.nam )
		self.int = INTERVAL( self.nam )
		
		self.gps_oss,self.sbas_oss,self.china_oss,self.gali_oss,self.obs_v2= TYPE_OBS( self ) 
		self.gps_week_ref, self.gps_sow_ref = GPSTIME( self.nam )
		self.data = ""
		self.phi  = "" ; self.lmbd = "" ; self.hell = ""
	#######################################################################  
	def COORD_GEOG(self):
		'''
		input:  self
		output: Latitude (F) and Longitude (L) and h WGS84
		'''
		# Longitudine calcolabile senza iterazioni
		x = self.xyz[0]; y = self.xyz[1]; z = self.xyz[2]
		L = np.arctan2(y,x)
		r = (x**2 + y**2)**0.5
		# WGS84 parameters
		a = 6378137
		f = 1/298.257223563
		b = a*(1-f)
		e = (1-(b**2)/(a**2))**0.5
		# I step
		# ip h=0
		h1 = 0.0
		F  = np.arctan2(z, r*(1-e**2))
		Rn = a/((1-(e**2)*(np.sin(F))**2))**0.5
		# Max 100 iterations --> usually break after 4 iterations
		for i in xrange(0,100):
				h  = r/np.cos(F) - Rn
				if abs(h - h1) <= 0.00001: break
				F  = np.arctan2((z*(Rn+h)),(Rn*(1-e**2)+h)*r)
				Rn = a/((1-(e**2)*(np.sin(F))**2))**0.5
				h1 = h
		L_grad = L/(np.pi)*180
		F_grad = F/(np.pi)*180

		self.phi  = F_grad #latitudine
		self.lmbd = L_grad #longitudine
		self.hell = h      #altezza ellissoidica
	########################################################

	def READ_RINEX(self):
		'''
		Function to stock the obs of the rinex
		'''
		import re
		import datetime
		L1 = 1.57542e9    # Hz
		L2 = 1.22760e9    # Hz
		L5 = 1.17645e9   #Hz E5a
		L5b = 1.207140e9 #Hz E5b
		L5a_b = 1.191795e9   #Hz E5a+b

		c  = 299792458    # m/s
		#
		lam1=c/L1         # m
		lam2=c/L2         # m
		lam5=c/L5         #m

		num_lin = 1

		if int(self.obs_v2[0])<=5:
				num_lin = 2	
		
		f   = open(self.nam, "r")
		lns = f.readlines()
		f.close()

		sats = []; ora  = []; sod  = []; c1   = []; l1   = []; l2   = []; s1=[];s2=[]; c5=[]; l5=[]
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
					#print ln
					if len(ln)==0: continue
					
					if len(ln[0])==2:
						#print ln
						if len(ln)>8:
							def tween(seq, sep): #The following will add a "separator" element between each of those in a list https://stackoverflow.com/questions/5920643/add-an-item-between-each-item-already-in-the-list
								return reduce(lambda r,v: r+[sep,v], seq[1:], seq[:1])
							ln[7:]=tween(ln[7:],'0') #add a zero
							ln[7:] = [reduce(lambda x, y: x + y, ln[7:])]
						else:
							ln=ln

						try:
							n_sats=int(ln[-1][0:2])
						except ValueError:
							n_sats=int(ln[-1][0])

						sats_appoggio1=[]
						sats_appoggio2=[]
						if n_sats<=9:
								for rr in np.arange(0,3*n_sats,3):
									satellite = ln[-1][1:][rr:rr+3]
									sats.append(satellite)

						else:
								if (n_sats/12.)%1!=0:#riga in piu'
									#satellite1 = ln[-1][2:][rr:rr+3]
									#sats_appoggio1.append(satellite1)
									riga=n_sats/12
									for tt in xrange(1,riga+1):
										for rr1 in np.arange(0,3*(12),3):
											if len(lns[i+tt].split())>1:
												ciao=lns[i+1].split()
												ciao=tween(ciao,'0')
												riga_da_convertire=[reduce(lambda x, y: x + y, ciao)]
												lns[i+tt]=" ".join(str(x) for x in riga_da_convertire)

												satellite2 = lns[i+tt].split()[0][rr1:rr1+3]
												sats_appoggio2.append(satellite2)
											else:
												satellite2 = lns[i+tt].split()[0][rr1:rr1+3]
												sats_appoggio2.append(satellite2)
									for rr2 in np.arange(0,3*n_sats,3):
										satellite1 = ln[-1][2:][rr2:rr2+3]
										sats_appoggio1.append(satellite1)

								else:
									riga=n_sats/12-1
									for tt in xrange(1,riga+1):
										for rr1 in np.arange(0,3*(12),3):
											if len(lns[i+tt].split())>1:
												ciao=lns[i+1].split()
												ciao=tween(ciao,'0')
												riga_da_convertire=[reduce(lambda x, y: x + y, ciao)]
												lns[i+tt]=" ".join(str(x) for x in riga_da_convertire)

												satellite2 = lns[i+tt].split()[0][rr1:rr1+3]
												sats_appoggio2.append(satellite2)
											else:
												satellite2 = lns[i+tt].split()[0][rr1:rr1+3]
												sats_appoggio2.append(satellite2)
									for rr2 in np.arange(0,3*n_sats,3):
										satellite1 = ln[-1][2:][rr2:rr2+3]
										sats_appoggio1.append(satellite1)

						satellite1=np.concatenate((np.asarray(sats_appoggio1),np.asarray(sats_appoggio2)))
						bau=np.where(map(lambda miao: len(miao)==0,satellite1))
						if len(bau[0])!=0:
								satellite1=np.delete(satellite1,bau[0])
						satellite=satellite1.tolist()
						for xw in xrange(0,len(satellite)):
									sats.append(satellite[xw])
								#map(lambda xw: sats.append(xw),satellite)
						#sats.append(satellite)


						if len(self.obs_v2[1:])<=5:

							for k in xrange(1, n_sats+1):
								if n_sats<=12:
									successiva=lns[i+k]
								else:
									sats_more=int(n_sats/12)
									successiva=lns[i+k+sats_more]
								dict_oss={}
								for g,h in zip(xrange(0,len(self.obs_v2[1:])),np.arange(1,81,16)):#196        										      
									try:
										dict_oss[str(self.obs_v2[1:][g])]=float(successiva[h:h+14])

									except ValueError:
										dict_oss[str(self.obs_v2[1:][g])]=float('Nan')

								l1_value_new  = (dict_oss['L1']) *lam1  
								l2_value_new = dict_oss.get('L2',float('Nan')) *lam2 
								c1_value_new = dict_oss.get('C1','P1')
								s1_value_new = dict_oss.get('S1',float('Nan'))
								s2_value_new = dict_oss.get('S2',float('Nan'))
								c5_value_new = dict_oss.get('C5',float('Nan'))
								l5_value_new = dict_oss.get('L5',float('Nan'))*lam5 
								
 								

								c1.append(   c1_value_new   )
								l1.append(   l1_value_new )
								l2.append(   l2_value_new  )
								c5.append(   c5_value_new  )
								l5.append(   l5_value_new  )

								s1.append(   s1_value_new  )
								s2.append(   s2_value_new  )

								timing = (  datetime.time( int(ln[3]), int(ln[4]) , int(float(ln[5][0:2])) ) )
								ora.append ( timing.isoformat() )
								sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))



						else:
							righez=int(len(self.obs_v2[1:]))/5 #righe in piu' rispetto prima riga
							if Decimal(len(self.obs_v2[1:])/5.)%1!=0:
								righe=int(len(self.obs_v2[1:]))/5								
							else:
								righe=int(len(self.obs_v2[1:]))/5-1
							for k in xrange(1, n_sats*(righe+1),righe+1):#for k in xrange(1, n_sats*2+1,2)
								#righe=int(len(self.obs_v2[1:]))/5
								dict_oss={}
								for qq, ee in zip(xrange(0,righe+1),np.arange(1,len(self.obs_v2[1:])+1,5)):
									if n_sats<=12:
										successiva=lns[i+k+qq]
										#print successiva
										#print i,k,qq,ee
									else:
										if Decimal(n_sats/12.)%1!=0:										    
											sats_more=int(n_sats/12)
										else:
											sats_more=int(n_sats/12)-1
										successiva=lns[i+k+qq+sats_more]
									for g,h in zip(xrange(0,len(self.obs_v2[ee:])),np.arange(1,81,16)):#196        										      
										try:
											dict_oss[str(self.obs_v2[ee:][g])]=float(successiva[h:h+14])
										except ValueError:
											dict_oss[str(self.obs_v2[ee:][g])]=float('Nan')


								
								l1_value_new  = (dict_oss['L1']) *lam1  
								l2_value_new =  (dict_oss['L2']) *lam2 
								c1_value_new = dict_oss.get('C1','P1')
								s1_value_new = dict_oss.get('S1',float('Nan'))
								s2_value_new = dict_oss.get('S2',float('Nan'))
								c5_value_new = dict_oss.get('C5',float('Nan'))
								l5_value_new = dict_oss.get('L5',float('Nan'))*lam5 
								
								

								c1.append(   c1_value_new   )
								l1.append(   l1_value_new )
								l2.append(   l2_value_new  )
								c5.append(   c5_value_new  )
								l5.append(   l5_value_new  )

								s1.append(   s1_value_new  )
								s2.append(   s2_value_new  )
							
								timing = (  datetime.time( int(ln[3]), int(ln[4]) , int(float(ln[5][0:2])) ) )
								ora.append ( timing.isoformat() )								
								sod.append(int(ln[3])*60*60 + int(ln[4])*60 + float(ln[5][0:3]))

		if len(sats)>len(sod):
			start_time = time.time()
			indici=[] # se ho dei numeri tra i satelliti, cosi li elimino
			for uu in xrange(0,len(sats)):
				try:
					if type(float(sats[uu]))==float:
						indici.append(uu)
				except ValueError:
						continue
			print "--- %s seconds for elimnate number from sat---" % (time.time() - start_time)

			if len(indici)!=0:
				sats=np.delete(sats,indici)
				sats=sats.tolist()
		else:
			pass


		sats_arr = np.asarray(sats)                                          
		ora_arr  = np.asarray(ora)
		sod_arr  = np.asarray(sod)
		c1_arr = np.asarray(c1)
		l1_arr = np.asarray(l1)
		l2_arr = np.asarray(l2)
		c2_arr = np.zeros(len(l1))

		c5_arr = np.asarray(c5)
		l5_arr = np.asarray(l5)

		s1_arr = np.asarray(s1)
		s2_arr = np.asarray(s2)
		self.data= sats_arr, ora_arr, sod_arr, l1_arr, l2_arr, c1_arr,c2_arr,s1_arr,s2_arr, c5_arr,l5_arr

