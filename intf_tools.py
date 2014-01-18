#!/usr/bin/python

"""intf_toos
The main shared library for interferometer processing type operations.  
Included toos for displaying, processing and simulating lightning 
interferometer type data.  
"""
########
# Written by Michael Stock, 2012
# versions: 
# the intf_tools version must match intf_process for processing to work!
#	5	-	Most of the individual functions broken out into a intf_tools library 
#			the global settings object introduced
#	6	-	Major revision to the file format, now including a header.  
#			Added capability for non perpendicular baselines, and more than 3 antennas
#	6.1	-	Moved the intf_process argument parser out of intf_tools (why was it there)
#			the Settings object less specific to intf_process

version = 6.1


#----------- Imports -------------
#numpy does all the math and calculations and supporting
from numpy import *
import numpy as np	#this covers the times when I do things the 'right' way
#os does file operations (chdir, etc) and can run external programs
#sys does program operations (exit, get argruments, etc)
#argparse parses command line arguments (once you get them)
import os,sys,argparse,struct,time,gzip
from matplotlib import colors


###################
# Globals

#the noise spectra, we initialize this as none, but a noise spectra 
#must be initialized for the facny windowing functions to work.
NoiseSpectra= { 	(0,1):None,
					(0,2):None,
					(0,3):None,
					(1,2):None,
					(1,3):None,
					(2,3):None,
					0: None,
					1: None,
					2: None,
					3: None }

#radius of the earth
Re   = 6371000		    #meters

#this is the global settings class
Settings   = None #to be initialized later

#~ colorDict = {	'red': ( 	(0.00, 0.3, 0.3),
							#~ (0.15, 0.0, 0.0),
							#~ (0.20, 0.0, 0.0),
							#~ (0.50, 0.0, 0.0),
							#~ (0.70, 1.0, 1.0),
							#~ (0.90, 0.9, 0.9),
							#~ (1.00, 1.0, 1.0) ),
				#~ 'green':(	(0.00, 0.0, 0.0),
							#~ (0.15, 0.0, 0.0),
							#~ (0.25, 0.3, 0.3),
							#~ (0.45, 1.0, 1.0),
							#~ (0.70, 0.9, 0.9),
							#~ (0.90, 0.0, 0.0),
							#~ (1.00, 0.0, 0.0) ),
				#~ 'blue':(	(0.00, 0.5, 0.5),
							#~ (0.20, 1.0, 1.0),
							#~ (0.45, 1.0, 1.0),
							#~ (0.50, 0.1, 0.1),
							#~ (0.70, 0.0, 0.0),
							#~ (0.90, 0.0, 0.0),
							#~ (1.00, 0.9, 0.9) ) }
#~ cmap_mjet   = colors.LinearSegmentedColormap('mjet',colorDict,256)

### same color map, goes to red at the end
#~ colorDict = {	'red': ( 	(0.00, 0.3, 0.3),
							#~ (0.15, 0.0, 0.0),
							#~ (0.20, 0.0, 0.0),
							#~ (0.50, 0.0, 0.0),
							#~ (0.70, 0.9, 0.9),
							#~ (1.00, 1.0, 1.0) ),
				#~ 'green':(	(0.00, 0.0, 0.0),
							#~ (0.15, 0.0, 0.0),
							#~ (0.25, 0.3, 0.3),
							#~ (0.45, 1.0, 1.0),
							#~ (0.70, 0.9, 0.9),
							#~ (1.00, 0.0, 0.0) ),
				#~ 'blue':(	(0.00, 0.5, 0.5),
							#~ (0.20, 1.0, 1.0),
							#~ (0.45, 1.0, 1.0),
							#~ (0.50, 0.1, 0.1),
							#~ (0.70, 0.0, 0.0),
							#~ (1.00, 0.0, 0.0) ) }
#~ cmap_mjet   = colors.LinearSegmentedColormap('mjet',colorDict,256)

### again same color map, goes to dark red at the end
colorDict = {	'red': ( 	(0.00, 0.3, 0.3),
							(0.15, 0.0, 0.0),
							(0.20, 0.0, 0.0),
							(0.50, 0.0, 0.0),
							(0.70, 1.0, 1.0),
							(0.90, 1.0, 1.0),
							(1.00, 0.7, 0.7) ),
				'green':(	(0.00, 0.0, 0.0),
							(0.15, 0.0, 0.0),
							(0.25, 0.3, 0.3),
							(0.45, 1.0, 1.0),
							(0.70, 0.9, 0.9),
							(0.90, 0.0, 0.0),
							(1.00, 0.0, 0.0) ),
				'blue':(	(0.00, 0.5, 0.5),
							(0.20, 1.0, 1.0),
							(0.45, 1.0, 1.0),
							(0.50, 0.1, 0.1),
							(0.70, 0.0, 0.0),
							(0.90, 0.0, 0.0),
							(1.00, 0.0, 0.0) ) }
cmap_mjet   = colors.LinearSegmentedColormap('mjet',colorDict,256)

colorDict = {	'red': ( 	(0.00, 0.0, 0.2),
							(0.30, 1.0, 1.0),
							(0.50, 1.0, 1.0),
							(0.60, 0.0, 0.0),
							(1.00, 0.0, 0.0) ),
				'green':(	(0.00, 0.0, 0.0),
							(0.30, 0.0, 0.0),
							(0.40, 0.3, 0.3),
							(0.50, 0.0, 0.0),
							(0.60, 0.3, 0.3),
							(1.00, 0.0, 0.0) ),
				'blue':(	(0.00, 0.0, 0.0),
							(0.40, 0.0, 0.0),
							(0.50, 1.0, 1.0),
							(0.70, 1.0, 1.0),
							(1.00, 0.0, 0.2) ) }
cmap_rb   = colors.LinearSegmentedColormap('rb',colorDict,256)

#####
# Classes
class SettingsClass:
	"""SettingsClass
	Store much of the information about processing and location of the 
	interferometer in a single object.
	A settings object is intended to exist as a global so that all 
	intf_tools functions have access to this information.
	"""
	
	def __init__(self):	
		"""SettingsClass.__init__(self)
		creates a settings object
		"""
		#this is the location of antenna 1, used for LMA comparisons, 
		#the measurement wasn't very acurate
		self.intfLoc   = [    33.97836 *pi/180,
							-107.19353 *pi/180,
							3163.5]
		
		#the antenna locations, in meters
		#these are based on tape measurements.  While they've been 
		#phased together based on 2013/07/08 data, they may still be off by 
		#a constant factor.
		self.antLocs    = array( [ [ -5.47,-15.39 ],
								   [  0.00,  0.00 ],
								   [-15.97,  5.86 ],
								   [-28.54,-15.63 ] ] )
		
		#Antennas delays (in ns).
		#phased together based on 2013/07/08 data
		self.antDels    = array( [ 	198.00,
									210.00,
									192.03,
									162.03 ] )
		########
		# Default Values (files)
		self.inFileS    = None
		self.outFileS   = None
		self.noiseFileS = None
		
		# Default Values (processing settings)
		self.numSamples = 256
		self.numIterate = 64
		self.numInterp  = 1
		self.windwoing  = None
		self.numAnts	= 3
		
		# Default Values (flags)		
		self.noise      = False	#all windows are assumed to be noise
		self.storeSpec  = False	#storing noise spectra information?
		self.debug 	    = False	#debugging output, not used
		self.hdf5       = False	#write an hdf5 file as well?

		# Default Values (Constants)
		# These values shouldn't change
		self.sampleRate = 180.		#the sample rate of the digitizer in MS/s
		self.bandwidth  = [20,80]
		self.c		    = 3.0e8		#speed of light
		self.minSep     = self.numSamples/2	#separation between 2 events, in samples
		
		# Defaults (Counters)
		self.SatCounter  = 0	#counts the occurrence of saturation
		self.startSample = 0	#were to start processing (0 is the beginning)
		self.stopSample  = -1	#were to stop processing (-1 is the end)

		# Defaults (Misc)
		self.aRange	    = None
		self.eRange     = None
		self.arguments  = None

	def file_initialize(self):
		"""SettingsClass.file_initialize(self)
		writes the header for the output file
		"""
		if not self.outFileS:
			#the path to the output file hasn't been set
			return
		#initialize the output file
		f = gzip.open(self.outFileS,'w')
		f.write( self.header.report() )
		f.write( self.report() )
		f.close()

	def file_finalize(self):
		"""SettingsClass.file_finalize(self)
		this goes through and updates some of the header values 
		which may have changed during processing"""
		if self.SatCounter == 0:
			#we have nothing to update
			return
		print 'finalizing file (this make take a second)'
		f = gzip.open(self.outFileS,'r')
		lines_i = f.readlines()
		f.close()
		for i in range(100):
			lineS = lines_i[i].strip()
			if lineS.split()[0] == '#Saturation':
				#update the number of saturations
				lines_i[i] = '#Saturation  '+ repr(self.SatCounter)		+ '\n'
				break
		f = gzip.open(self.outFileS,'w')
		f.writelines(lines_i)
		f.close()	
		
		#write an hdf5 copy of the data?
		if self.hdf5:
			#read in the data as an object
			data = read_data_file(self.outFileS)
			if self.outFileS[-6:] == 'dat.gz':
				hd5FileS = self.outFileS[:-6]+'hdf5'
			else:
				hd5FileS = self.outFileS[:-4]+'.hdf5'
			#write this to an hdf5 file
			write_hdf5_file(data, hd5FileS)
	
	def report(self):
		"""SettingsClass.report()
		prints all the settings"""
		
		output = ""
		
		output += '#inFileS     :'+ repr(self.inFileS)			+ '\n'
		output += '#outFileS    :'+ repr(self.outFileS)			+ '\n'
		output += '#noiseFileS  :'+ repr(self.noiseFileS)		+ '\n'
		output += '#startSample :'+ repr(self.startSample)		+ '\n'
		output += '#verbose     :'+ repr(self.verbose)			+ '\n'
		output += '#debug       :'+ repr(self.debug)			+ '\n'
		output += '#windowing   :'+ repr(self.windowing)		+ '\n'
		output += '#noise       :'+ repr(self.noise)			+ '\n'
		output += '#storeSpec   :'+ repr(self.storeSpec)		+ '\n'
		output += '#sampleRate  :'+ repr(self.sampleRate)		+ '\n'
		output += '#bandwidth   :'+ repr(self.bandwidth)		+ '\n'
		output += '#numSamples  :'+ repr(self.numSamples)		+ '\n'
		output += '#numIterate  :'+ repr(self.numIterate)		+ '\n'
		output += '#numInterp   :'+ repr(self.numInterp)		+ '\n'
		output += '#c           :'+ repr(self.c)				+ '\n'
		output += '#############'+ '\n'
		output += '#AntsUsed    :'+ repr(self.numAnts)				 + '\n'
		output += '#Ant0loc     :'+ repr(list(self.antLocs[0]))  + '\n'
		output += '#Ant0delay   :'+ repr(self.antDels[0])        + '\n'
		output += '#Ant1loc     :'+ repr(list(self.antLocs[1]))  + '\n'
		output += '#Ant1delay   :'+ repr(self.antDels[1])        + '\n'
		output += '#Ant2loc     :'+ repr(list(self.antLocs[2]))  + '\n'
		output += '#Ant2delay   :'+ repr(self.antDels[2])        + '\n'
		output += '#Ant3loc     :'+ repr(list(self.antLocs[3]))  + '\n'
		output += '#Ant3delay   :'+ repr(self.antDels[3])        + '\n'
		output += '#############'+ '\n'
		output += '#minSep      :'+ repr(self.minSep)			 + '\n'
		output += '#Saturation  :'+ repr(self.SatCounter)		 + '\n'
		#output += '#arguments  ', self.arguments
		output += '#############'+ '\n'
		output += '#Time, Azimuth, Elevation, cosa, cosb, Pk2Pk, RMS, eXpk, eCls, eStd, eMlt, Red, Green, Blue, startSample, iMax \n'
		return output
		
class RawHeader:
	"""RawHeader
	This stores and parses the header information in the raw data files
	"""

	def __init__(self,fileS):
		
		####
		# Default Values
		#the most important things the header has
		#version info
		version		= 0
		size		= 0
		#time info
		trigTime	= 0
		year        = 0
		month 		= 0
		day 		= 0
		hour 		= 0
		minute 		= 0
		second 		= 0
		usecond     = 0
		TrigType    = 0
		
		#data info
		sampleRate  = 0
		chAEnabled  = False
		chAACCoupled= False
		chAImpedence= False
		chABWLimit  = False
		chARange    = 0
		chBEnabled  = False
		chBACCoupled= False
		chBImpedence= False
		chBBWLimit  = False
		chBRange    = 0
		
		version, size, header = self.read(fileS)
		
		#this is the most important thing
		self.version 	= version
		self.size 		= size
		
		#now we parse the header data
		if version == 1 or version == 2:
			(self.trigTime,
			 self.sampleRate,
			 self.decimation,
			 self.autoDmaFlags,
			 self.samplesPerRecord,
			 self.recordsPerTransfer,
			 self.transfersPerAcquisition,
			 self.preTriggerSamples,
			 self.channelMask,
			 self.chARange,
			 self.chBRange) = struct.unpack('<Q7LBLL',header)
			#I don't care if preTriggerSamples is set, in this 
			#version pretrigger didn't work
			self.preTriggerSamples = 0
			#decode the time information
			date = time.gmtime(self.trigTime)
			self.year   = date[0]
			self.month  = date[1]
			self.day    = date[2]
			self.hour   = date[3]
			self.minute = date[4]
			self.second = date[5]
			#decode the channel mask
			channelMask = self.channelMask
			if (channelMask >> 0) % 2:
				self.chAEnabled = True
			else:
				self.chAEnabled = True
			if (channelMask >> 1) % 2:
				self.chAACCoupled = True
			else:
				self.chAACCoupled = False
			if (channelMask >> 2) % 2:
				self.chAImpedence = 50
			else:
				self.chAImpedence = 1e6
			if (channelMask >> 3) % 2:
				self.chABWLimit = True
			else:
				self.chABWLimit = False
			if (channelMask >> 4) % 2:
				self.chBEnabled = True
			else:
				self.chBEnabled = False
			if (channelMask >> 5) % 2:
				self.chBACCoupled = True
			else:
				self.chBACCoupled = False
			if (channelMask >> 6) % 2:
				self.chBImpedence = 50
			else:
				self.chBImpedence = 1e6
			if (channelMask >> 7) % 2:
				self.chBBWLimit = True
			else:
				self.chBBWLimit = False

		elif version >= 5 and version <= 12:
			#headVersion
			#HeadSize
			#year
			#Month
			#Day
			#hour
			#minute
			#second
			#uSecond
			#TrigType
			#SampleRate
			#Decimation
			#Samples/Record
			#ChRange
			#CouplingId
			#ImpedanceId
			#BWLimit
			(self.year,
			 self.month,
			 self.day,
			 self.hour,
			 self.minute,
			 self.second,
			 self.usecond,
			 self.TrigType,
			 self.sampleRate,
			 self.decimation,
			 self.samplesPerRecord) = struct.unpack('<6Hd4I16x',header)
			self.recordsPerTransfer = None
			self.transfersPerAcquisition = None
			if self.TrigType == 0:
				#External Trigger
				self.preTriggerSamples = int(0.5*self.sampleRate)
			elif self.TrigType == 1:
				#Manual Trigger #1
				self.preTriggerSamples = 5*self.sampleRate
			elif self.TrigType == 2:
				#Manual Trigger #1
				self.preTriggerSamples = 10*self.sampleRate
			else:
				self.preTriggerSamples = 0
		else:
			print version
			raise VersionError

	def read(self,fileS):
		"""RawHeader.read(fileS)
		reads the raw header, but doesn't parse anything other than the 
		version and size
		"""
		#test data file and returns the header if it is valid, 
		#returns None otherwise
		f = open(fileS,'r')
		(version,size) = struct.unpack('HH',f.read(4))
		head = f.read(size-4)
		f.close()
		return version, size, head

	def report(self):
		"""RawHeader.report()
		prints all the settings"""

		output = ""
		
		TrigTimeS = "%04d/%02d/%02d %02d:%02d:%02d.%06i"%(
					self.year,self.month,self.day,
					self.hour,self.minute,self.second, self.usecond)
		
		output += '#DataFileVersion    :'+ repr(self.version) 	 + '\n'
		output += '#TriggerTime        :'+ repr(TrigTimeS) 		 + '\n'
		output += '#Year               :'+ repr(self.year)       + '\n'
		output += '#Month              :'+ repr(self.month)      + '\n'
		output += '#Day                :'+ repr(self.day)        + '\n'
		output += '#Hour               :'+ repr(self.hour)       + '\n'
		output += '#Minute             :'+ repr(self.minute)     + '\n'
		output += '#Second             :'+ repr(self.second)     + '\n'
		output += '#uSecond            :'+ repr(self.usecond)    + '\n'
		output += '#############'+ '\n'
		output += '#TriggerType        :'+ repr(self.TrigType)   + '\n'
		output += '#SampleRate         :'+ repr(self.sampleRate) + '\n'
		output += '#PreTriggerSamples  :'+ repr(self.preTriggerSamples) + '\n'
		output += '#############'+ '\n'
		
		return output
	
class Event:
	"""Event
	The Event class holds information about a single interferometer 
	source, it does not process the source however
	"""
	def __init__(self,startSample,eamp0,eamp1,expeak,eclos,cosa,cosb,ants,iMax,(r,g,b)):
		
		Az, El = cosab2AzEl_perp(cosa, cosb)
		#time
		self.startSample = startSample
		self.ants   = ants
		self.time   = (iMax-Settings.preTriggerSamples)/1000./Settings.sampleRate	#in ms
		self.iMax   = iMax
		self.sTime  = ( '%4.4f'%self.time).rjust( 9 )
		self.color  = (	str(r).rjust(3),
						str(g).rjust(3),
						str(b).rjust(3))
		self.ants   = ants	#the antenna indexes
		#then we store an array of past solutions
		self.iXpeak = [expeak]
		self.iAmp0  = [eamp0]
		self.iAmp1  = [eamp1]
		self.iClos  = [eclos]
		w = 1-(eclos**4)/(eclos**4+2.0**4)	#the weight
		self.cosa   = [[cosa,w]]
		self.cosb   = [[cosb,w]]

		#deal with the permutations
		self.permutations = 1.
		if Settings.numAnts == 4:
			self.permutations = 4.

		#the strings
		eclos = abs(eclos)
		self.strings( Az,El,cosa,cosb,expeak,eamp0,eamp1,eclos,100,1 )
		
			
	def printout(self):
		print "%s %s %s %s %s %s %s %s %s %s %s %s %s %s  %i  %i"%\
			(	self.sTime,
				self.Az,
				self.El,
				self.sCosa,
				self.sCosb,
				self.sAmp0,
				self.sAmp1,
				self.sXpeak,
				self.sclos,
				self.sstd,
				self.mult,
				self.color[0],
				self.color[1],
				self.color[2],
				self.startSample, 
				self.iMax-self.startSample)
	
	def write(self,fileS):
		f = gzip.open(fileS,'a')
		f.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s  %i  %i\n"%\
			(	self.sTime,
				self.Az,
				self.El,
				self.sCosa,
				self.sCosb,
				self.sAmp0,
				self.sAmp1,
				self.sXpeak,
				self.sclos,
				self.sstd,
				self.mult,
				self.color[0],
				self.color[1],
				self.color[2],
				self.startSample, 
				self.iMax-self.startSample) )
		f.close()
	def strings(self,Az,El,cosa,cosb,expeak,eamp0,eamp1,eclos,estd,mult ):
		#calculate the fractional multiplicity
		mult /= float(Settings.numSamples)/Settings.numIterate*self.permutations
		#Errors
		self.sCosa  = ('%1.4f'%cosa     ).rjust( 8 )
		self.sCosb  = ('%1.4f'%cosb     ).rjust( 8 )
		self.sXpeak	= ('%0.2f'%expeak	).rjust( 5 )
		self.sAmp0	= ('%4.1f'%eamp0	).rjust( 6 )
		self.sAmp1	= ('%4.1f'%eamp1	).rjust( 6 )
		self.sstd	= ('%3.2f'%estd		).rjust( 7 )
		self.sclos  = ('%2.2f'%eclos	).rjust( 6 )
		self.mult   = ('%0.2f'%mult     ).rjust( 5 )
		#location
		self.Az  	= ('%3.3f'%Az		).rjust( 8 )
		self.El  	= ('%2.3f'%El		).rjust( 7 )
		
	def append(self,eamp0,eamp1,expeak,eclos,cosa,cosb):
		#append the erros
		self.iXpeak.append(expeak)
		self.iAmp0.append( eamp0 )
		self.iAmp1.append( eamp1 )
		self.iClos.append( eclos )
		#append the locations
		w = 1-(eclos**4)/(eclos**4+2.0**4)	#the weight
		self.cosa.append( [cosa,w] )
		self.cosb.append( [cosb,w] )
		#print array(self.iAz)
		
		#calculate Az and El
		solsa = array(self.cosa)
		solsb = array(self.cosb)
		
		cosa = sum( solsa[:,0]*solsa[:,1] )/sum(solsa[:,1])
		cosb = sum( solsb[:,0]*solsb[:,1] )/sum(solsb[:,1])
		
		#calculate the variance
		vara = sum( solsa[:,1]*(solsa[:,0]-cosa)**2 )/sum(solsa[:,1])
		varb = sum( solsb[:,1]*(solsb[:,0]-cosb)**2 )/sum(solsb[:,1])
		estd = sqrt( vara + varb )
		if estd < 1:
			estd = arcsin( estd )*180/pi
		else:
			estd = 100
		
		#calculate azimuth and elevation
		Az, El = cosab2AzEl_perp( cosa, cosb )
		
		#calculate the closure
		#note, the weight (solsa[:,1]) is based on the closure, so 
		#we are wieghting closure based on closure, which is funny
		eclos = sum( abs(array(self.iClos))*solsa[:,1] )/sum(solsa[:,1])

		#calculate new strings
		self.strings( Az,El,cosa,cosb,max(self.iXpeak),
					  max(self.iAmp0),max(self.iAmp1),eclos,
					  estd,len(self.cosa) )

class ProcHeader():
	
	def __init__(self,dic):
		
		#the dictionary can be generated with parse_ascii_header
		#or it can be pulled from the hdf5 file (I hope)
		self.dic = dic
		
		#see if we have a processed file version in here
		if 'ProcessedFileVersion' not in dic.keys():
			#that means it's v6.0
			dic['ProcessedFileVersion'] = 6.0
				
		for key in dic.keys():
			if isinstance( dic[key], np.ndarray ):
				dic[key] = dic[key].tolist()
			setattr(self,key,dic[key])

class LmaData():
	def __init__(self,fileS):
		self.fileS = fileS
		if fileS[-2:] == 'gz':
			f = gzip.open(fileS,'r')
		else:
			f = open(fileS,'r')
			
		#read forward to the data
		lineS = f.readline()
		while lineS.strip() != '*** data ***':
			lineS = f.readline()
		
		self.azRange = [0,360]	#degrees
		self.elRange = [0,90]
		self.caRange = [-1.,1.]	#cosine projection
		self.cbRange = [-1.,1.]
		self.xRange  = [-20,20]	#km
		self.yRange	 = [-20,20]
		self.zRange  = [2.8,20]
		self.ciRange = [0,100]	#chisq for lma data
		
		self._time = []
		self._azim = []
		self._elev = []
		self._rang = []
		self._lat  = []
		self._lon  = []
		self._x    = []
		self._y    = []
		self._z    = []
		self._chisq= []
		self._pwr  = []
		self._charge=[]
		self._mask = []
		self._cosa = []
		self._cosb = []

		lineS = f.readline().strip().split()
		t0 = int(float(lineS[0]))
		#t0 = 0
		self.t0 = t0	#unlikely to be used
		while lineS != []:
			if len(lineS) == 7:
				t    = float(lineS[0])
				lat  = float(lineS[1])*pi/180
				lon  = float(lineS[2])*pi/180
				alt  = float(lineS[3])
				chisq= float(lineS[4])
				pwr  = float(lineS[5])
				mask = int(lineS[6],16)
				charge=0
				stats= 0#not really, but I'm not ready do deal with this yet
			elif len(lineS) == 8:
				t    = float(lineS[0])
				lat  = float(lineS[1])*pi/180
				lon  = float(lineS[2])*pi/180
				alt  = float(lineS[3])
				chisq= float(lineS[4])
				stats= int(  lineS[5])
				pwr  = float(lineS[6])
				charge=0
				mask = int(  lineS[7], 16)				
			elif len(lineS) == 9:
				t    = float(lineS[0])
				lat  = float(lineS[1])*pi/180
				lon  = float(lineS[2])*pi/180
				alt  = float(lineS[3])
				chisq= float(lineS[4])
				stats= int(  lineS[5])
				pwr  = float(lineS[6])
				charge=int(  lineS[7])
				mask = int(  lineS[8], 16)
			else:
				print len(lineS)
			
			(Az,El,r) = calc_range_bearing(lat,lon,alt)
			t += r/3e8
			
			
			#calculate x,y,z
			(x,y,z) = calc_xyz( Az, El, r )
			
			self._time.append(t)
			self._azim.append(Az)
			self._elev.append(El)
			self._rang.append(r)
			self._lat.append(lat*180/pi)
			self._lon.append(lon*180/pi)
			self._x.append(x)
			self._y.append(y)
			self._z.append(z)
			self._chisq.append(chisq)
			self._pwr.append(pwr)
			self._charge.append(charge)
			self._mask.append(mask)
			self._cosb.append( cos(El/180*pi )*sin(Az/180*pi) )
			self._cosa.append( cos(El/180*pi )*cos(Az/180*pi) )

			lineS = f.readline().strip().split()

		f.close()
		self._time = array(self._time)
		self._azim = array(self._azim)
		self._elev = array(self._elev)
		self._rang = array(self._rang)
		self._lat  = array(self._lat)
		self._lon  = array(self._lon)
		self._alt  = array(self._z)
		self._x    = array(self._x)
		self._y    = array(self._y)
		self._z    = array(self._z)
		self._chisq= array(self._chisq)
		self._charge=array(self._charge)
		self._pwr  = array(self._pwr)
		self._mask = array(self._mask)
		self._cosa = array(self._cosa)
		self._cosb = array(self._cosb)
		
	
		#convert the time array
		self._time -= t0
		self._time *= 1000

		self.tRange = self._time[0],self._time[-1]
		print self.tRange

		#print self._time[0]
		self.time = []
		self.limits()
		
	def limits (self):
		""" limit the data based on the ranges"""
		mask = 	(self._time>=self.tRange[0] )&(self._time<=self.tRange[1]) &\
				(self._azim>=self.azRange[0])&(self._azim<=self.azRange[1]) &\
				(self._elev>=self.elRange[0])&(self._elev<=self.elRange[1]) &\
				(self._cosa>=self.caRange[0])&(self._cosa<=self.caRange[1]) &\
				(self._cosb>=self.cbRange[0])&(self._cosb<=self.cbRange[1]) &\
				(self._chisq>=self.ciRange[0])&(self._chisq<=self.ciRange[1])
		
		print len(self.time),
		self.time = self._time[mask]
		self.azim = self._azim[mask]
		self.elev = self._elev[mask]
		self.rang = self._rang[mask]
		self.lat  = self._lat [mask]
		self.lon  = self._lon [mask]
		self.alt  = self._alt [mask]
		self.x    = self._x   [mask]
		self.y    = self._y   [mask]
		self.z    = self._z   [mask]
		self.chisq= self._chisq[mask]
		self.charge= self._charge[mask]
		self.pwr  = self._pwr [mask]
		self.cosa = self._cosa[mask]
		self.cosb = self._cosb[mask]
		print len(self.time)
				
		

class ProcData():
	
	def __init__(self, header,data):
		
		########
		# The 2 main arrays
		self.mask    = None #this won't be None for long
		self.data    = None #also won't be None for long
		
		#we store the raw data (un filtered)
		#and we store the filtered data
		self.rawData = data
		if isinstance(header,str):
			#this is old old v3 data!!!
			try:
				#fake some important header items
				year   = int( header.split('_')[1].split('.')[0] )
				month  = int( header.split('_')[1].split('.')[1] )
				day    = int( header.split('_')[1].split('.')[2] )
				hour   = int( header.split('_')[2].split('.')[0] )
				minute = int( header.split('_')[2].split('.')[1] )
				second = int( header.split('_')[2].split('.')[2] )
			except:
				raise ValueError('could not parse filename: %s'%header)
				return
			header = {}
			header[ 'year'   ] = year
			header[ 'month'  ] = month
			header[ 'day'    ] = day
			header[ 'hour'   ] = hour
			header[ 'minute' ] = minute
			header[ 'second' ] = second
			header['TriggerTime'] = '%04i/%02i/%02i %02i:%02i:%02i'%(
									year,
									month,
									day,
									hour,
									minute,
									second)
			header['columns'] = [	'Time',
									'Azimuth',
									'Elevation',
									'cosa',
									'cosb',
									'Pk2Pk',
									'RMS',
									'eXpk',
									'eCls',
									'fwhm',
									'eStd',
									'eMlt',
									'Red',
									'Green',
									'Blue',
									'startSample',
									'iMax']
			#these are guesses
			header[ 'numIterate'] = 64
			header[ 'numSamples'] = 256
			header[ 'numInterp' ] = 2
			v3 = True
		if isinstance(header,dict):
			print 'building header from dictionary'
			header = ProcHeader(header)
		self.header = header
		
		########
		# These are some paraments that get used a lot
		self.nbins 	  = 200
		self.TriggerTimeS = header.TriggerTime
		self.tStart   = np.floor( min( self.rawData[0] ) )
		self.tStop    = np.ceil(  max( self.rawData[0] ) )
		self.tOffset  = 0
		#limits
		self.tRange   = [self.tStart,self.tStop]	#time
		self.azRange  = [0,360]		#azimuth
		self.elRange  = [-30,90]	#elevation
		self.caRange  = [-1.1,1.1]	#cosa
		self.cbRange  = [-1.1,1.1]	#cosb
		#quality
		self.tCls     = 1.25
		self.tStd     = 1.0
		self.tXpk     = 0.7
		self.tMlt     = 0.7

		#find the index of various values
		#if the header is proper, all these will be there
		self.iTime = header.columns.index('Time')
		self.iAzim = header.columns.index('Azimuth')
		self.iElev = header.columns.index('Elevation')
		self.icosa = header.columns.index('cosa')
		self.icosb = header.columns.index('cosb')
		self.ipkpk = header.columns.index('Pk2Pk')
		self.irms  = header.columns.index('RMS')
		self.ieXpk = header.columns.index('eXpk')
		self.ieCls = header.columns.index('eCls')
		self.ieStd = header.columns.index('eStd')
		self.ieMlt = header.columns.index('eMlt')
		self.iRed  = header.columns.index('Red')
		self.iBlue = header.columns.index('Green')
		self.iGrn  = header.columns.index('Blue')
		self.isSmp = header.columns.index('startSample')
		self.iiMax = header.columns.index('iMax')
			
		
		if self.header.numIterate < self.header.numSamples:
			#we do a normal filter
			print 'Overlapping Windows, normal Filter'
			self.filter()
		else:
			#there's no overlap, normal filtering won't work
			print 'Non-overlapping Windows, coarse Filter'
			self.tStd = 200
			self.tXpk = 0.6
			self.tMlt = max( self.rawData[self.ieMlt] ) +1
			self.filter()
		
		#self.update()
		#initialize the mask
		
		#make a hist of the rawData
		self.rawDataHist = np.histogram2d( self.elev,
										   self.time, 
										   weights=self.pkpk, 
										   bins=[self.nbins/2,self.nbins], 
										   range=([0,90],[self.tStart,self.tStop]) )		
	def time_from_trigger(self):
		"""time_from_trigger(self)
		adjusts the time array to count from the trigger time.  
		This is useful is you are trying to identify a period of time 
		to process, since intf_process counts from the trigger
		"""
		tOffset = 0
		self.tRange[0] = self.tRange[0] - self.tOffset + tOffset
		self.tRange[1] = self.tRange[1] - self.tOffset + tOffset
		self.rawDataHist[2][:] += tOffset - self.tOffset		
		self.tOffset   = tOffset
		self.update()
		return self.tOffset	#this is incase other programs code needs to track this
		
	
	def time_from_second(self):
		"""time_from_second(self)
		adjusts the time array to count from the beginning of the second, 
		or at least attempts to do that.  It is accurate to about 20 ms
		"""
		print self.tStart
		try:
			tStartOff = self.tStart-self.tStart%1000
			print 'tStart', self.tStart, tStartOff
			tOffset   = self.header.uSecond/1000-tStartOff
			print 'found usecond data', self.header.uSecond
		except:
			tOffset = -self.tStart
			print 'didn\'t find usecond data'
		if self.tRange[0] - self.tOffset + tOffset > 1000:
			tOffset -= 1000
		self.tRange[0] = self.tRange[0] - self.tOffset + tOffset
		self.tRange[1] = self.tRange[1] - self.tOffset + tOffset
		print 'tb',self.rawDataHist[2][0]
		self.rawDataHist[2][:] += tOffset - self.tOffset
		print 'tb',self.rawDataHist[2][0]
		
		self.tOffset   = tOffset
		self.update()
		return self.tOffset	#this is incase other programs code needs to track this
		
	
	def sort(self, arr):
		self.mask = self.mask[ arr.argsort() ]
		self.update()
	
	def reset_limits(self):
		self.tStart   = np.floor( min( self.rawData[0] ) )
		self.tStop    = np.ceil(  max( self.rawData[0] ) )
		self.tRange   = [self.tStart+self.tOffset,self.tStop+self.tOffset]	#time
		self.azRange  = [0,360]		#azimuth
		self.elRange  = [0,90]	#elevation
		self.caRange  = [-1.1,1.1]	#cosa
		self.cbRange  = [-1.1,1.1]	#cosb
		self.limits()
	
	def limits(self):
		#construct a new mask based on the data limts
		#the mask is a set of indices so that the data can be sorted 
		#in interesting ways
		self.mask = where(
					(self.data[self.iTime]>=self.tRange[0] - self.tOffset)&
					(self.data[self.iTime]<=self.tRange[1] - self.tOffset)&
					(self.data[self.iAzim]>=self.azRange[0])&
					(self.data[self.iAzim]<=self.azRange[1])&
					(self.data[self.iElev]>=self.elRange[0])&
					(self.data[self.iElev]<=self.elRange[1])&
					(self.data[self.icosa]>=self.caRange[0])&
					(self.data[self.icosa]<=self.caRange[1])&
					(self.data[self.icosb]>=self.cbRange[0])&
					(self.data[self.icosb]<=self.cbRange[1]) )[0]
		self.sort( self.mask )


	def update(self):
		#set all data variables
		#Time, Azimuth, Elevation, cosa, cosb, Pk2Pk, RMS, eXpk, eCls, eStd, eMlt, Red, Green, Blue, startSample, iMax   
		self.time = self.data[self.iTime][self.mask] + self.tOffset
		self.azim = self.data[self.iAzim][self.mask]
		self.elev = self.data[self.iElev][self.mask]
		self.cosa = self.data[self.icosa][self.mask]
		self.cosb = self.data[self.icosb][self.mask]
		self.pkpk = self.data[self.ipkpk][self.mask]
		self.rms  = self.data[self.irms ][self.mask]
		self.eXpk = self.data[self.ieXpk][self.mask]
		self.eCls = self.data[self.ieCls][self.mask]
		self.eStd = self.data[self.ieStd][self.mask]
		self.eMlt = self.data[self.ieMlt][self.mask]
		self.freq = self.data[self.iRed:self.iBlue,self.mask]/255. #this is stored as a color
		self.sSam = self.data[self.isSmp][self.mask]
		self.iMax = self.data[self.iiMax][self.mask]
		
		self.aMax = self.data[self.ipkpk].max()
		N = len(self.data[self.ipkpk])
		self.a95  = self.data[self.ipkpk][self.data[self.ipkpk].argsort()[ 99*N/100 ]]
		self.a05  = self.data[self.ipkpk][self.data[self.ipkpk].argsort()[  1*N/20  ]]
		self.aMin = self.data[self.ipkpk].min()
		print N,self.aMin, self.a05, self.a95, self.aMax
		
		#make the bghist for the overview
		self.dataHist = np.histogram2d(   self.elev, 
										  self.time, 
										  weights=self.pkpk, 
										  bins=[self.nbins/2,self.nbins], 
										  range=([0,90],[self.tStart,self.tStop]) )
					
	def filter( self ):
		#Store the threshold information
		eCls = self.tCls
		eStd = self.tStd
		eXpk = self.tXpk
		eMlt = self.tMlt
		
		mask         = ((self.rawData[self.ieCls]/eCls)**2+(self.rawData[self.ieStd]/eStd)**2 <  1) &\
					   ((self.rawData[self.ieXpk]/eXpk)**2+(self.rawData[self.ieMlt]/eMlt)**2 >  1)
		#mask = mask & (self.rawData[self.ieMlt] < eMlt)
		self.data = self.rawData[:,mask]

		self.limits()
		self.update()
		print 'loaded %i of %i values between %f and %f'%(
				len(self.time),len(self.rawData[0]),self.tStart,self.tStop)


	def coarse_filter(self, pkpk=None):
		if pkpk==None:
			pkpk = np.median( self.rawData[self.ipkpk] )*1.5
		mask = (self.rawData[self.ipkpk]>pkpk)
		self.data = self.rawData[:,mask]
		self.limits()
		self.update()


#########
# Data Functions


def read_raw_file_data(fileS = None, header=None, startSample = None, numSamples=None):
	"""read_raw_file_data(fileS = None, header=None, startSample = None, numSamples=None):
	Reads in the 4 channels of raw data, 
	includes handling for different raw data versions
	"""
	#get the variables we need in the local namespace
	if fileS == None:
		fileS       = Settings.inFileS
	if header == None:
		header      = Settings.header
	if startSample == None:
		startSample = Settings.startSample
	if numSamples == None:
		numSamples  = Settings.numSamples
	#for version 2, the data is interwoven in blocks
	if header.version == 2:
		#the data is stored in block of 1024*1024 samples on 2 channels
		#just in case, read the number of samples per block from the header, 
		#multiply by 4 (2 channels, 2 words per sample)
		blockSize = header.samplesPerRecord*header.recordsPerTransfer*4
		startLoc  = startSample*4
		#if startLoc is bigger than a block, add 1 blockSize for ever
		#blocksize bigger it is
		startLoc += int(startLoc/blockSize)*blockSize
		#print startSample,startSample*4,blockSize,startLoc
		f = open(fileS)
		#seek to start location
		try:
			f.seek(startLoc+header.size)
			chAB = struct.unpack('HH'*numSamples,f.read(numSamples*4))
			f.seek(startLoc+blockSize+header.size)
			chCD = struct.unpack('HH'*numSamples,f.read(numSamples*4))
			f.close()
		except:
			#the read failed for some reason, probably EOF
			return None
		#now shove the data into an array
		data = array([
			array(chAB).reshape(numSamples,2)[:,0],
			array(chAB).reshape(numSamples,2)[:,1],
			array(chCD).reshape(numSamples,2)[:,0],
			array(chCD).reshape(numSamples,2)[:,1]])
	#the alternative is Mark's data, which is split into 4 files
	elif header.version >= 5 and header.version <= 12:
		try:
			startLoc  = int(startSample*2)
			#channel A
			f = open( fileS[:-3]+'chA','r')
			f.seek(startLoc+header.size)
			chA = struct.unpack( '%iH'%numSamples, f.read(numSamples*2) )
			f.close()
			#channel B
			f = open( fileS[:-3]+'chB','r')
			f.seek(startLoc+header.size)
			chB = struct.unpack( '%iH'%numSamples, f.read(numSamples*2) )
			f.close()
			#channel C
			f = open( fileS[:-3]+'chC','r')
			f.seek(startLoc+header.size)
			chC = struct.unpack( '%iH'%numSamples, f.read(numSamples*2) )
			f.close()
			#channel D
			f = open( fileS[:-3]+'chD','r')
			f.seek(startLoc+header.size)
			chD = struct.unpack( '%iH'%numSamples, f.read(numSamples*2) )
			f.close()
		except:
			return None

		data = array( [
			array( chA ),
			array( chB ),
			array( chC ),
			array( chD ) ])
	return data

def read_raw_waveform_file_data( fileS = None, header=None, startSample=None, numSamples=None, maxSamples=50000.):
	"""read_raw_waveform_file_data( fileS = None, header=None, startSample=None, numSamples=None, maxSamples=10000):
	reads in data from a single channel for use as plotting a waveform
	
	Only works for version 5 and above data
	
	maxSamples - 
		If more than maxSamples data is requested, the waveform 
		is decimated
	"""
	blockSize  = 512	#how many samples to read at a time
	maxSamples = float(maxSamples)	#needs to be a float for later math
	#get the variables we need in the local namespace
	if fileS == None:
		fileS       = Settings.inFileS
	if header == None:
		header      = Settings.header
	if startSample == None:
		startSample = Settings.startSample
	if numSamples == None:
		numSamples  = Settings.numSamples

	#determine the decimation amount
	decimation = 0
	while decimation == 0:
		if numSamples > maxSamples:
			decimation = int( maxSamples*blockSize/numSamples )
		else:
			decimation = blockSize
		if decimation == 0:
			blockSize *=2

	wave = []
	time = []
	#version 5 and above data, split into files
	if header.version >= 5 and header.version <= 12:
		iSample = startSample
		#how many samples/block should we keep
		iTime = arange(blockSize)
		if not os.path.exists(fileS):
			return None
		while iSample < startSample + numSamples:		
			try:
				startLoc  = int(iSample*2)
				#channel whatever we're on
				f = open( fileS,'r')
				f.seek(startLoc+header.size)
				chunk = struct.unpack( '%iH'%blockSize, f.read(blockSize*2) )
				f.close()
			except:
				#hit end of file
				break
			chunk = array(chunk)
			#calculate the decimation mask
			mask  = (abs(chunk-mean(chunk))).argsort()[-decimation:]
			mask.sort()
			#append to the previous data
			wave  += chunk[mask].tolist()
			time  += (iTime[mask]+iSample).tolist()
			
			#increment counter
			iSample += blockSize
	else:
		return None
	#truncate
	output = array( [time,wave],dtype='float' )
	print output.shape, decimation
	output = output[ :, output[0,:]<= startSample+numSamples ]
	#convert time to ms
	#(iMax-Settings.preTriggerSamples)/1000./Settings.sampleRate
	output[0,:] = (output[0,:]-header.preTriggerSamples)*1000./header.sampleRate
	return  output

def write_hdf5_file(data, outFileS):
	import h5py
	#does the file exit?
	if os.path.exists(outFileS):
		print 'removing output file'
		#remove it
		os.remove(outFileS)
	#header might be a header object, or it might just be the dictionary
	#we'd rather just work with the object
	with h5py.File(outFileS, 'w') as f:
		#create the dataset (the big things)
		dset = f.create_dataset("data", data.rawData.shape, dtype='f')
		dset[...] = data.rawData
		#set the attributes
		for key in sorted( data.header.dic.keys() ):
			if data.header.dic[key] == None:
				dset.attrs[key] = 'None'
			else:
				dset.attrs[key] = data.header.dic[key]

def read_data_file(inFileS):
	"""Determines filetype and reads"""
	extS = os.path.splitext(inFileS)[1]
	if (extS == '.hdf5') or (extS == '.h5'):
		data = read_hdf5_file(inFileS)
	elif (extS == '.dat') or (inFileS[-6:]=='dat.gz'):
		data = read_ascii_file(inFileS)
	elif extS == '.csv':
		data = read_ascii_file(inFileS)
	else:
		raise TypeError('Invalid File Type')
	
	return data
		
def read_hdf5_file(inFileS):
	import h5py
	with h5py.File(inFileS,'r') as f:
		dset = f['data']
		header,data = parse_hdf5_data(dset)
	return ProcData(header,data)
	
def read_ascii_file(inFileS):
	if inFileS[-2:] == 'gz':
		f = gzip.open(inFileS)
	else:
		f = open(inFileS)
	#read in all the lines
	lines_i = f.readlines()
	f.close()
	
	header,data = parse_ascii_data(lines_i)
	if header != None:
		return ProcData(header,data)
	else:
		return ProcData(inFileS,data)
		
def parse_ascii_header(lines_i):
		output = {}
		for lineS in lines_i:
			if not ':' in lineS:
				#this is a spacer
				continue
			i = lineS.index(':')
			labelS = lineS[1:i].strip()
			valueS = lineS[i+1:].strip()
			#labelS is the key (easy)
			#the value though needs to be converted
			if valueS[0] == "'":
				#this is a string, just strip those characters
				valueS = valueS[1:-1]
			#next are the flags
			elif valueS == 'None':
				valueS = None
			elif valueS == 'True':
				valueS = True
			elif valueS == 'False':
				valueS = False
			elif valueS[0] == '[':
				#a list, this is a bit more complex
				valueS = valueS[1:-1].split(',')
				for j in range(len(valueS)):
					valueS[j] = float(valueS[j].strip())
			#all that's left are ints and floats
			elif '.' in valueS:
				#it's a float
				valueS = float(valueS)
			else:
				#it's an int
				valueS = int(valueS)
			#put the value into the dic
			output[labelS] = valueS
		
		#the last line should be the column values
		columns = lines_i[-1][1:].split(',')
		for i in range(len(columns)):
			columns[i] = columns[i].strip()
		output['columns'] = columns
		
		return output

def parse_ascii_data(lines_i):
	#build the header
	i = 0
	while lines_i[i][0] == '#':
		i+=1
	if i > 0:
		header = ProcHeader( parse_ascii_header(lines_i[:i]) )
	else:
		header = None

	#build the data array (the big thing)
	#initialize the array
	H = i
	N = len(lines_i)-i
	M = len(lines_i[i].split())
	#M better be 16
	print N,M
	if M < 16:
		#old v3 data without cosa, cosb
		M += 2
	data = np.empty( (M,N) )
	#go through each line and convert
	while i < N+H:
		lineS = lines_i[i]
		if header == None:
			#v3 data
			tmp  = np.fromstring(lineS, sep=' ')
			cosa = cos(tmp[1]*pi/180)*cos(tmp[2]*pi/180)
			cosb = sin(tmp[1]*pi/180)*cos(tmp[2]*pi/180)
			data[:3,i-H] = tmp[:3]
			data[3 ,i-H] = cosa
			data[4 ,i-H] = cosb
			data[5:,i-H] = tmp[3:]
		else:
			data[ :,i-H ] = np.fromstring(lineS, sep=' ')
		i += 1
	
	return header, data

def parse_hdf5_data(dset):
	header = ProcHeader( dict( dset.attrs) )
	data = dset[...]
	
	return header,data
	
def read_noise_spectra(fileS=None):
	global NoiseSpectra
	if fileS == None:
		fileS = Settings.noiseFileS
	
	if not os.path.exists(fileS):
		print "Error: noise spectra missing"
		sys.exit(1)
	spec = loadtxt(fileS)
	NoiseSpectra[ (0,1) ] = spec[:,0] + spec[:,1]*1J
	NoiseSpectra[ (0,2) ] = spec[:,2] + spec[:,3]*1J
	NoiseSpectra[ (0,3) ] = spec[:,4] + spec[:,5]*1J
	NoiseSpectra[ (1,2) ] = spec[:,6] + spec[:,7]*1J
	NoiseSpectra[ (1,3) ] = spec[:,8] + spec[:,9]*1J
	NoiseSpectra[ (2,3) ] = spec[:,10]+ spec[:,11]*1J
	NoiseSpectra[   0   ] = spec[:,12]
	NoiseSpectra[   1   ] = spec[:,13]
	NoiseSpectra[   2   ] = spec[:,14]
	NoiseSpectra[   3   ] = spec[:,15]

def init_noise_spectra(fileS=None):
	global NoiseSpectra

	for key in NoiseSpectra.keys():
		NoiseSpectra[key] = zeros( Settings.numSamples )
	
def write_noise_spectra(fileS=None):
	if fileS == None:
		fileS = Settings.noiseFileS
	spec = zeros([2*Settings.numSamples,16])
	spec[:,0] = real( NoiseSpectra[ (0,1) ][0] )
	spec[:,1] = imag( NoiseSpectra[ (0,1) ][0] )
	spec[:,2] = real( NoiseSpectra[ (0,2) ][0] )
	spec[:,3] = imag( NoiseSpectra[ (0,2) ][0] )
	spec[:,4] = real( NoiseSpectra[ (0,3) ][0] )
	spec[:,5] = imag( NoiseSpectra[ (0,3) ][0] )
	spec[:,6] = real( NoiseSpectra[ (1,2) ][0] )
	spec[:,7] = imag( NoiseSpectra[ (1,2) ][0] )		
	spec[:,8] = real( NoiseSpectra[ (1,3) ][0] )
	spec[:,9] = imag( NoiseSpectra[ (1,3) ][0] )		
	spec[:,10]= real( NoiseSpectra[ (2,3) ][0] )
	spec[:,11]= imag( NoiseSpectra[ (2,3) ][0] )		

	spec[:,12]= real( NoiseSpectra[   0   ][0] )
	spec[:,13]= real( NoiseSpectra[   1   ][0] )
	spec[:,14]= real( NoiseSpectra[   2   ][0] )
	spec[:,15]= real( NoiseSpectra[   3   ][0] )
	
	print 'Saving noise spectra averaged over %i windows to:'%NoiseSpectra[ 0 ][1], fileS
	savetxt( fileS,spec)	

#####
# Functions

def pkpk(arr):
	return arr.max() - arr.min()

def rms(arr):
	return (mean(arr**2))**.5

def calc_range_bearing(lat,lon,alt):
	#distance
	#d is along the earth
	a = sin( (Settings.intfLoc[0]-lat)/2 )**2 + cos(Settings.intfLoc[0])*cos(lat)*sin( (Settings.intfLoc[1]-lon)/2 )**2
	d = (Re)*2*arctan2( sqrt(a), sqrt(1-a) )

	#print d, alt2, alt1
	#elevation
	El = arctan( abs(alt-Settings.intfLoc[2])/d )
	#get the range
	r  = sqrt( d**2 + abs(alt-Settings.intfLoc[2])**2 )

	#bearing
	Az = arctan2( sin(lon-Settings.intfLoc[1])*cos(lat), cos(Settings.intfLoc[0])*sin(lat)-sin(Settings.intfLoc[0])*cos(lat)*cos(lon-Settings.intfLoc[1]) )
	if Az < 0:
		Az += 2*pi
	
	return (Az*180/pi, El*180/pi, r)

def calc_xyz( Az, El, r ):
	z = r*sin( El/180*pi ) + Settings.intfLoc[2]
	d = r*cos( El/180*pi )
	x = d*sin( Az/180*pi )
	y = d*cos( Az/180*pi )
	return x,y,z


def cosab2AzEl(cosa, cosb, ants):
	theta = antTheta( ants )
	basel = antBaseline(ants)
	
	dtheta = theta[0]-theta[1]
	
	Az = theta[0] + arctan2( 
		(cosa * cos(dtheta) - cosb),
		(cosa * sin(dtheta))
		)
	
	Az = rad2deg( Az )
	
	cosEl = sqrt( 
		(cosa**2 + cosb**2 - 2*cosa*cosb*cos(dtheta)) / sin(dtheta)**2 
		)
	if cosEl <= 1:
		El = rad2deg( arccos( cosEl ) )
	elif cosEl <= 2:
		#it's probably just beyond the horizon, which might be ok
		#use negative angles to handle this
		El = -rad2deg( arccos( 2-cosEl ) )
	else:
		#how'd this get here?!
		El = 0
	
	return Az, El
def cosab2AzEl_perp(cosa, cosb):
	"""same as cosab2AzEl but it assumes perpendicular, non-rotated 
	frame"""
	
	Az = arctan2( cosb, cosa )
	Az = rad2deg( Az )
	
	cosEl = sqrt( cosa**2 + cosb**2)
	if cosEl <= 1:
		El = rad2deg( arccos( cosEl ) )
	elif cosEl <= 2:
		#it's probably just beyond the horizon, which might be ok
		#use negative angles to handle this
		El = -rad2deg( arccos( 2-cosEl ) )
	else:
		#how'd this get here?!
		El = 0
	
	return Az, El

def cosab2cosab(cosa, cosb, ants):
	"""convert direction cosines from non-perpendicular rotated baselines 
	to direction cosines from perpendicular non-rotated baselines"""
	
	#cosa = cos(Az)cos(El)
	#cosb = sin(Az)cos(El)
	
	#find Azimuth	
	theta = antTheta( ants )
	basel = antBaseline(ants)
	
	dtheta = theta[0]-theta[1]
	
	Az = theta[0] + arctan2( 
		(cosa * cos(dtheta) - cosb),
		(cosa * sin(dtheta))
		)
	
	#find cos(El)
	cosEl = sqrt( 
		(cosa**2 + cosb**2 - 2*cosa*cosb*cos(dtheta)) / sin(dtheta)**2 
		)
	
	#recalculate cosa and cosb for perpendicular baselines
	#note, cosEl might be larger than 1 here
	cosa  = cosEl*cos(Az)
	cosb  = cosEl*sin(Az)
	
	return cosa, cosb

def antTheta( ants ):
	output = []
	for i in range(1, len(ants)):
		dLoc = Settings.antLocs[ants[0]]-Settings.antLocs[ants[i]]
		output.append( 
			arctan2( dLoc[0], dLoc[1] )
			)
	return output

def antBaseline( ants ):
	output = []
	for i in range(1, len(ants)):
		output.append( 
			sqrt( sum( 
				(Settings.antLocs[ants[0]] - Settings.antLocs[ants[i]])**2 
			)  )
		)
	return output

def rad2deg( rad ):
	"""converts radians to degrees
	requires the conversion to be between 0-360 (positive angle)"""
	rad *= 180/pi
	if rad < 0:
		rad += 360
	return rad

def impulse(t, w, N):
	"""def impulse(t, w, N):
	t	- time of the pulse (in samples)
	w	- width of the pulse (in samples)
	N	- length of the window (in samples)
	---
	out - array"""
	t = float(t)
	iN = arange(N)
	envelope = exp( -.5*((iN-t)/w)**2 )
	out_i = random.randn( N )*envelope
	
	return out_i

def filter( arr, f_nyq, bw, df=3 ):
	"""def filter( arr, f_nyq, bw ):
	arr		- the array to filter
	f_nyq	- the nyquist frequency
	bw		- array, the frequency to filter to
	---
	out		- array"""
	#df = 3	#this is how long it takes to go from low to high
	f_min = 0	#minimum of the filter
	f0 = [ bw[0]-df, bw[0]+df ]
	f1 = [ bw[1]-df, bw[1]+df ]
	N = arr.shape[0]
	
	freq = abs( fft.fftfreq(2*N,1./f_nyq) )
	#build the filter
	W = zeros( (2*N, arr.shape[1]) )
	
	#the rise
	mask =  (freq>=f0[0])&(freq<=f0[1])
	for i in range(arr.shape[1]):
		W[mask,i] = 0.5*(1-cos( pi*(freq[mask]-f0[0])/(f0[1]-f0[0]) ))
	#the fall
	mask =  (freq>=f1[0])&(freq<=f1[1])
	for i in range(arr.shape[1]):
		W[mask,i] = 0.5*(1-cos( pi*(freq[mask]-f1[1])/(f1[1]-f1[0]) ))
	#the top
	W[ (freq>f0[1])&(freq<f1[0]),: ] = 1
	#the bottom
	W[ W<f_min ] = f_min
	
	#apply the fitler
	out_i = fft.fft( arr, 2*N, axis=0 )
	out_i *= W
	out_i = fft.ifft( out_i, axis=0 )
	
	return real( out_i[:N] )

def ddelay( arr, t ):
	"""def delay( arr, t ):
	This funciton seems to have run into a bug in python
	Large negative values of t cause issues
	"""
	out = arr.copy()
	t = -t
	if t > 0:
		out[:-t]  = out[t:]
		out[-t:] = 0
	if t < 0:
		out[-t:]  = out[:t]
		out[:-t]  = 0
	
	return out

def fdelay( arr, t):
	"""fdelay( arr, t):
	delay arr by (float) t samples"""
	
	N = len(arr)
	
	#fft, pad by a factor of 2 to help with circular issues
	out = fft.fft(arr, 2*N)
	
	#calculate the delay factor
	df = exp( -2J*pi*fft.fftfreq(2*N)*t )
	out *= df
	out = fft.ifft(out)[:N]
	
	return real(out)
	
def decimate( arr, N):
	"""decimate( arr, N):
	decimate arr to length N"""
	dN = len(arr)/float(N)
	mask = round_( arange( N ) * dN ).astype(int)
	mask[ mask> len(arr)-1 ] = len(arr)-1
	out_i = arr[mask]
	
	return out_i

def xcorr( d0, d1, P):
	"""xcorr( d0, d1, P)
	Computes the cross correlation of d0 and d1
	"""
	N = len(d1)
	
	fft0 = fft.fft(d0, 2*N)
	fft1 = fft.fft(d1, 2*N)
	
	#multiply
	x01 = fft0 * (fft1.conj())
	#interpolate
	x01 = fpad(x01,P)
	#invert
	x01 = fft.fftshift( fft.ifft(x01) )[1:]
	#normalize
	x01 /= sqrt( fft.ifft(fft0*fft0.conj())[0]*fft.ifft(fft1*fft1.conj())[0] )

	return real(x01)

def fwhm( x ):
	"""fwhm(x)
	Calculates the full width half max of array x in a quick and dirty 
	sort of way.  Most importantly it's quick
	"""
	x = abs( x-mean(x) )
	xMax = max(x)
	iX = where( x>xMax/2 )[0]
	return iX[-1]-iX[0]

def fpad(X,M):
	"""fpad(X,M)
	Frequency 0 pads X to be M*len(X) long,
	Used for fft based interpolation
	input - 
		X	-	fft of a signal, X(f)
		M	-	factor to interpolate by
	output -
		padded output(f) 
	"""
	
	if M  <= 1:
		return X
	N = len(X)
	output = append( X[:N/2], X[N/2]/2. )
	output = append( output, zeros((M-1)*N-1) )
	output = append( output, X[N/2]/2. )
	output = append( output, X[N/2+1:] )
	output *= M	#the scaling is needed so amplitude of the ifft remains unchanged
	
	return output
	
def xmax(x):
	xMax = argmax(x)
	#these catches should never come up
	if xMax == 0:
		xMax +=1
	if xMax == len(x)-1:
		xMax -=1
	para_fit = polyfit([0,1,2],x[xMax-1:xMax+2],2)
	#take the derivative
	para_max = [para_fit[0]*2,para_fit[1]]
	para_max = roots(para_max)
	#the happens if you put all 0's into the max finder
	if not para_max:
		return 0
	else:
		return xMax + para_max[0]

def store_spectra(ffts,ants):
	"""store_spectra(ffts, ants)
	appends the data in ffts to the noise spectra
	each spectra in NoiseSpectra has the average spectra information and 
	number of windows which contributed to make the average.  This way 
	the average can be updated
	input - 
		ffts	-	a tuple of 3 fft's
		ants	-	a tuple of the antenna numbers
	output - 
		None
	"""
	global NoiseSpectra
	
	#go through the permutations of the cross spectra
	i = 0
	while i < len(ants) -1:
		j = i+1
		while j < len(ants):
			key = tuple( sorted( (ants[i], ants[j]) ) )
			
			if NoiseSpectra[ key ] == None:
				#then this spectra hasn't been set yet
				NoiseSpectra[key] = (ffts[i]*ffts[j].conj(),1)
			else:
				#it has been set, so we average this fft with the last one
				N = float(NoiseSpectra[key][1])
				NoiseSpectra[key] = (NoiseSpectra[key][0]*(N/(N+1)) + \
										ffts[i]*ffts[j].conj()/(N+1),
										N+1)
			j += 1
		i += 1
	#store the auto spectra
	for i in range(len(ants)):
		key = ants[i]
		if NoiseSpectra[key] == None:
			#then this spectra hasn't been set yet
			NoiseSpectra[key] = ( abs(ffts[i])**2, 1)
		else:
			#it has been set, so we average this fft with the last one
			NoiseSpectra[key] = (NoiseSpectra[key][0]*(N/(N+1)) + \
									abs(ffts[i])**2/(N+1),
									N+1)

def symm_smooth( arr,M ):
	"""symm_smooth( arr, M)
	Symmetrically smooth array arr with and order M filter
	Used for estimating spectra for the windowing functions
	(this function could use improvement)
	
	input - 
		arr	-	array
		M	-	order
	ouput - 
		array
	"""
	#smooth in a symmetric way and returns an array of the same length
	for i in range(M):
		arr = (arr + append( 0,arr[:-1]) + append( arr[1:],0 ))/3
	return arr

def calc_delay(d1,d2):
	"""calc_delay
	This is a simplified version of the process_window function 
	in intf_process.  It does the delay calculation on any array 
	2 arrays.  This is useful both for simulations, and for testing 
	new algorythms
	
	input - 
		d1	-	the first array
		d2	-	the second array
	output -
		
	"""
	####
	# STEP 1 Preliminary Stuff
	N = len(d1)
	if Settings != None:
		I     = gSettings.numIterate
		P     = gSettings.numInterp
		bw    = gSettings.bandwidth
		sRate = gSettings.SampleRate
	else:
		I     = N
		P     = 4
		bw    = [20,80]
		sRate = 180

	#weight the data towards the center of the window
	Wt = sin( arange( N ) *pi /N )*2./3+1./3
	d1 = ( d1 - mean(d1) )*Wt
	d2 = ( d2 - mean(d2) )*Wt
	
	#calculate the ffts
	fft1 = fft.fft(d1,2*N)
	fft2 = fft.fft(d2,2*N)
	
	freq = abs( sRate*fft.fftfreq( 2*N ) )
	W = zeros(2*N)+1
	#the bandpass filter is super super simple
	W[ freq<bw[0] ] = 0	#The stuff outside the hardware filters
	W[ freq>bw[1] ] = 0
	W0 = W.copy()	#we save a copy of this fitler if we use fancy filtering

	####
	# STEP 2 Metrics 1
	#there are some metrics that we need to calculate on d1 and d2
	#before we calculate the cross correlation
	
	
	#apply the filter to the data
	d1 = real( fft.ifft( fft1*W0 )[:N] )
	d2 = real( fft.ifft( fft2*W0 )[:N] )	
	
	#iMax is used for event differentiation, it is the location 
	#of the maximum event in the arrays
	iMax = [argmax(abs(d1)),
			argmax(abs(d2)) ]
	

	####
	# STEP 3 Cross Correlation
	
	#Multiply
	x12 = W * fft1 * (fft2.conj())
	#Interpolate
	x12 = fpad(x12,P)
	#Invert
	x12 = fft.fftshift( fft.ifft(x12) )[1:]
	#normalize (divide by the auto-corr at 0)
	x12 /= fft.ifft( fft1*fft1.conj() )[0] * fft.ifft( fft2*fft2.conj() )[0]
	#Drop Imaginary Part
	x12 = real( x12 )
	
	#calculate the time shift
	s12 = (xmax(x12)/P -N)/sRate
	
	####
	# STEP 4 Metrics 2
	epkpk = min([ d1.max() - d1.min(), d2.max() - d2.min() ])
	erms  = min([ rms(d1), rms(d2) ])
	expk  = max( x12 )
	
	return s12, iMax, [epkpk,erms,expk]

def place_label(ax,s,ar=1,size=12,box=True,figSize=(6.5,4)):
	import matplotlib as mpl
	#the ar is used because for axes with equal aspect, the 
	#get_position call returns bogus shit
	offset = 0.07	#in inches!
	texth  = size/100.
	#this is more complicated than it should be
	axBox = ax.get_position().get_points()
	axAsp = ax.get_data_ratio()
	w = diff(axBox[:,0])[0]*figSize[0]
	h = diff(axBox[:,1])[0]*figSize[1]*ar
	print 'w,h',w,h
	bbox=mpl.patches.Rectangle( (0,1-(2*offset+texth)/h),
			width=(2*offset+texth)/w,
			height=(2*offset+texth)/h,
			facecolor=(1,1,1),
			alpha=.8,
			transform=ax.transAxes)
	#an odd way of doing this
	if box:
		p = ax.add_patch(bbox)
		p.set_zorder=10
	t = ax.text( 0+offset/w,1-(offset+texth)/h,s , 
		horizontalalignment='left',
		verticalalignment='baseline',
		transform=ax.transAxes,
		size=size)
	t.set_zorder(11)


##########
# Initialize some things
Settings = SettingsClass()
