#!/usr/bin/python
# \file VARION.py
#  This script implements the VARION algorithm for real-time 
#  detection of sTEC variations using GNSS observations.
# author Giorgio Savastano, 23.10.2015. giorgio.savastano[AT]gmail.com
# modified by Michela Ravanelli, michela.ravanelli[AT]uniroma1.it
#
# Notice: Please acknowledge the use of the above software in any publications:
#    ``VARION software was provided by G. Savastano et al.,
#      and is available at URL: https://github.com/giorgiosavastano/VARION''.
#
# Reference: Savastano, G. et al. Real-Time Detection of Tsunami Ionospheric Disturbances with a Stand-Alone
# GNSS Receiver: A Preliminary Feasibility Demonstration. Sci. Rep. 7, 46607; doi: 10.1038/srep46607 (2017).
#
# Please send a copy of such publications to either G. Savastano or A. Komjathy:
#  Giorgio Savastano                      		     	Dr. Attila Komjathy
#  Civil, Building and Environmental Engineering,    	Ionospheric and Atmospheric Remote Sensing Group, 
#  University of Rome "La Sapienza"           		 	Jet Propulsion Laboratory, California Institute of Technology,
#  Rome, Italy.                               		 	Pasadena, California, USA.
#  E-mail: giorgio.savastano[AT]uniroma1.com         	E-mail: attila.komjathy[AT]jpl.nasa.gov 
# 
# -------------------------------------------------------------------------
#
# Copyright (C) 2015-2020  (see AUTHORS file for a list of contributors)
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
#------------------------------------------------------------------------
## IMPORT MODULES AND CLASSES ##
import argparse  
import os                           
import glob
import numpy as np
#
import myObs as mO
import myFunc as mF
import mySatsFunc as mSF
import RinexClass as RC

import pdb

parser = argparse.ArgumentParser(prog="VARION.py", description= """VARION.py is a script that process RINEX obs files
									   and apply the VARION algorithm in order to obtain sTEC measuraments.
									  author: Giorgio Savastano - giorgiosavastano@gmail.com 
 									  author:  Michela Ravanelli - michela.ravanelli@uniroma1.it """)

parser.add_argument("-staz", type=str, nargs='*', default="all", dest="stazName", help="This argument determines the station(s) will be processed." \
									  " By default, this parameter is set to process all the RINEX observation files in the working folder (/obs). ")
						   
parser.add_argument("-time", nargs='*', type=str, default="all", dest="analysisTime", 
										help="If no argument is given, the analysis is executed for " \
									  " all the time vector of the obs file." \
									  " Otherwise, the argument refers to the time for which the analysis"\
									  " should be performed and has to be in the format hh:min (GPS time)"\
									  "(e.g., 18:34 19:00)")
									  
parser.add_argument("-sat", type=str, nargs='*', default=0, dest="satNumber", help="This argument determines the satellite(s) will be considered." \
									  "By default, this parameter is set to process all the satellites in view for each epochs."\
									  "(e.g., G01 G05 G23)")    

parser.add_argument('-brdc', dest="brdcOrb",  action='store_true', help="This argument set the processing with the brdc file")       

parser.add_argument('-height', type=int, default=300, dest="hIono",  help='This argument determines the ionospheric shell height'\
										'By default, this value is set to 300 km')             
									   
print "VARION is a free and open source software that processes RINEX obs files in order to estimate sTEC values."
print "author: Giorgio Savastano - giorgiosavastano@gmail.com "
print "author: Michela Ravanelli - michela.ravanelli@uniroma1.it"
print " "
########################################################
## CLASSES ##
class myStation:
	""" This is the simple class to describe a GNSS station """
	def __init__ (self):
		self.name = ""
		self.oFile = ""        ## RINEX obs files
		self.GPSnFile = ""     ## RINEX GPS nav file
		self.GLOnFile = ""     ## RINEX GLO nav file
		self.skyFile  = ""     ## VADASE sky file
		self.brdcFile = ""
		self.process_able = False
	def VADASE_PROCESS_ABLE(self):
		""" This function checks if the station has the observation and sky
			file available in the processing folder. If this is the 
			case the "process_able" variable is set to True """
		if os.path.isfile(self.oFile) and os.path.isfile(self.GPSnFile):
			self.process_able = True
		elif os.path.isfile(self.brdcFile):
			self.process_able = True
########################################################
## VARIABLES ##                                       
#### Constant ####
L1 = 1.57542e9                           # HZ
L2 = 1.22760e9                           # HZ
A  = 40.308e16
c = 299792458.0                          # m/s

const_tec = ((L1**2)*(L2**2))/(A*(L1**2-L2**2))
sats = np.asarray( ['G01','G02','G03','G04','G05','G06','G07','G08','G09','G10','G11','G12',\
					   'G13','G14','G15','G16','G17','G18','G19','G20','G21','G22','G23','G24',\
					   'G25','G26','G27','G28','G29','G30','G31'] )
######################################################## 
os.chdir('..')
main_dir = os.getcwd()
obs_dir  = main_dir + '/obs'
out_dir  = main_dir + '/outputs'
os.chdir('obs')

## PROGRAM STARTS ##
args = parser.parse_args()
print args 
h_iono = args.hIono * 1000.0             # height of the ionospheric layer

if args.stazName == "all":
	stations = glob.glob('*.??o')
	stations.sort()
else:
	statio = args.stazName
	suffix = glob.glob(statio[0] + '*.??o')[0][4:]

	stations = [ sta + suffix for sta in statio ]
	stations.sort()
##########################################################   
#### METTERE POI OPZIONE PER IL FILE BRDC E IGS       ------- > IMPORTANTE 
if args.brdcOrb == True:
	brdc_file = glob.glob ( 'brdc' + '*.??n')
	print brdc_file

## COUNT HOW MANY NAVIGATION FILES ARE NOT AVAILABLE ##
myStationsProc = []                                                     # List of stations that will be processed

for sFile in stations:
		alreadyThere = False
		for station in myStationsProc:
			if sFile[0:4] == station.name:
				## The station is already in the list                   #
				## check if it has the observation and sky file         #
				## and if not assign them                               #
				if args.brdcOrb == True:
					station.brdcFile = brdc_file[0]
				if not station.oFile:
					station.oFile = sFile
				if not station.GPSnFile:
					sGPSnFile = sFile[:-1] + 'n' 
					if os.path.isfile(sGPSnFile):
						station.GPSnFile = sGPSnFile
				
				station.VADASE_PROCESS_ABLE()
							
				alreadyThere = True
				break
		## The station is not in the list
		
		if not alreadyThere:
			sStation = myStation()
			sStation.name = sFile[0:4]
			sStation.oFile = sFile
			sGPSnFile = sFile[:-1] + 'n' 
			if os.path.isfile(sGPSnFile):
				sStation.GPSnFile = sGPSnFile
			if args.brdcOrb == True:
				sStation.brdcFile = brdc_file[0]
			
			sStation.VADASE_PROCESS_ABLE()	
			myStationsProc.append(sStation)	

for i in myStationsProc:
		print i.name, i.oFile, i.GPSnFile, i.brdcFile, i.process_able 
##########################################################
if args.analysisTime != "all":
	 start = int(args.analysisTime[0][:2])*60.0*60.0 + int(args.analysisTime[0][3:5])*60.0   
	 stop  = int(args.analysisTime[1][:2])*60.0*60.0 + int(args.analysisTime[1][3:5])*60.0
	 print start, stop
	 
if args.satNumber == 0:
	sats_write = sats
	print sats_write
else:
	sats_write = np.asarray(args.satNumber)
	sats_write.sort()
	print sats_write

################################################################################
## EXECUTE VARION ##	

info_file = open(  "info.txt" , "w" )
for i in myStationsProc:
	if i.process_able:

		if args.brdcOrb == True:
			rinex_nav = brdc_file[0]
		else:
			rinex_nav = i.GPSnFile
		# CREATE THE RINEX OBJECT FROM THE CLASS RinexFile()
		rinex_obs = RC.RinexFile( i.oFile )

		lat_g,lon_g, h = mF.coord_geog(  rinex_obs.xyz[0],rinex_obs.xyz[1],rinex_obs.xyz[2]  )
		
		info_file.write( str(rinex_obs.nam)+ "\t" + str(rinex_obs.int) + "\t" + str(lat_g) + "\t" + str(lon_g) + "\n"  )
		##
		# read the rinex with the method built inside the class
		import time
		start_time = time.time()
		rinex_obs.READ_RINEX()
		rinex_obs.COORD_GEOG()
		print "RINEX file %s has been read in" % rinex_obs.nam
		print("--- %s seconds ---" % (time.time() - start_time))

		rinex_obs.data  =    mSF.skip_nan(rinex_obs,rinex_obs.data[5])

		#select just the satellite in view
		sats_write_1 = mSF.sat_selection( rinex_obs, sats_write, start, stop ) 
		
		try:
			start_time = time.time()
			sIP = mSF.coord_satellite( rinex_nav, rinex_obs, sats_write_1)
			print "Coord satellites has been computed in"
			print("--- %s seconds ---" % (time.time() - start_time))
			
		except ValueError:
			print 'station ' + str(rinex_obs.nam) + ' has been skipped'
			continue
################################################################################
		lista_G = []
		sIP_G_list = []
		data_list = []
		start_time = time.time()	
		for sa in sats_write_1:
				varion = mO.obs_sat( rinex_obs.data[0], rinex_obs.data[1], rinex_obs.data[2], rinex_obs.data[3], rinex_obs.data[4], sa )
				data_list.append( rinex_obs.data )  
				lista_G.append( varion )
				sIP_sat = mSF.track_sat( sIP, sa, start, stop  )
				####
				phi_ipp, lambda_ipp, h_ipp = mSF.coord_ipps( rinex_obs.xyz[0],rinex_obs.xyz[1],rinex_obs.xyz[2], sIP_sat[2], sIP_sat[3], sIP_sat[4], h_iono)
				sIP_G_list.append(  (sIP_sat[0],sIP_sat[1],phi_ipp,lambda_ipp,sIP_sat[6],sIP_sat[7])  )
		print "VARION algorithm has been computed and"
		print "IPP location has been computed for the satellites selected in"
		print("--- %s seconds ---" % (time.time() - start_time))

		################################################################################
		### REMOVE THE OUTLAYER
		stec_list = []
		sod_list = []
		for i in xrange( len(sats_write_1) ):
				mask = mF.no_outlayer_mask( lista_G[i][0] * const_tec / rinex_obs.int )  ## modify the treshold to remove the outlayer
				stec_list.append(  lista_G[i][0][mask] * const_tec / rinex_obs.int  ) 
				sod_list.append(  lista_G[i][2][mask]  )
		  
		################################################################################
		### POLINOMIAL INTERPOLATION OF THE DATA
		X_list = []; Y_list = []
		mask_list = []
		diff_list = []
		cum_list  = []
		import warnings
		warnings.simplefilter('ignore', np.RankWarning)
		for i in xrange( len(sats_write_1) ):
				X = sod_list[i]
				Y = stec_list[i] 
				mask       = (X>=start) & (X<=stop)
				try:
					p        = np.poly1d(  np.polyfit(X[mask], Y[mask], 10)  )
					interpo = p(  X[mask]  )
					# residual
					diff = Y[mask] - interpo  
					# integrate
					cum = mF.integrate(  diff, rinex_obs.int  )
					# append
					X_list.append(X)
					Y_list.append(Y)
					mask_list.append(mask)
					diff_list.append(diff)
					cum_list.append(cum)
				except (TypeError, IndexError):
					X_list.append(0.0)
					Y_list.append(0.0)
					mask_list.append(0.0)
					diff_list.append(0.0)
					cum_list.append(0.0)
		################################################################################
		### Create the .txt file
		################################################################################
		for i in xrange( len(sats_write_1) ):
				mask = (sIP_G_list[i][0] >= start) & (sIP_G_list[i][0] <= stop)
		
				f = open(out_dir + '/' + str( rinex_obs.nam ) +'_' + str(sats_write_1[i]) + '_' + str(args.hIono) + '.txt', 'w')
				f.write('sod' + '\t' + '\t'  + '\t' + 'sTEC' + '\t' + '\t'+ '\t' 'lon' + '\t' + '\t'+ '\t' 'lat'+ '\t' + '\t'+ '\t' 'ele' + '\t' + '\t'+ '\t' 'azi' '\n')
				try:
					for k in xrange( len(cum_list[i]) ):
						try:
							#### FIX DIFF OF TIME BETWEEN COORDINATES AND STEC (ONE COME FROM NAVIGATION FILE THE OTHER FROM OBS)
							## BUG FIXED  --> try with 30 s data
							inde = (np.where(X_list[i][mask_list[i]][k] ==  sIP_G_list[i][0][mask]) )
							f.write( str(sIP_G_list[i][0][mask][inde[0][0]]) + '\t' + '\t' + str(cum_list[i][k]) + '\t' + '\t' + \
										str(sIP_G_list[i][3][mask][inde[0][0]]) + '\t' + '\t' + str(sIP_G_list[i][2][mask][inde[0][0]]) + \
										    '\t' + '\t' + str(sIP_G_list[i][-1][mask][inde[0][0]])+'\t' + '\t' +str(sIP_G_list[i][4][mask][inde[0][0]]) +'\n')
						except IndexError:
							continue
				except TypeError or IndexError:
					continue
				f.close() 
info_file.close()  
	