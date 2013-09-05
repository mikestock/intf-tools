#!/usr/bin/python
#######
# intf_process.py <FILE_IN> [FILE_OUT]
#
#	Reads FILE_IN and calculates azimuth/elevation for the data
#	Input is AlazarTech 9642 4 channel digitized data, should be 
#	version 2.
#	If FILE_OUT is supplied, writes ascii azimuth elevation to 
#	that file, otherwise writes to FILE_IN.csv
#
# Written by Michael Stock, 2012
# versions:
#	1	-	Initial code, sort of worked
#	3	-	Robust code, windowing functions implemented
#	4	-	Attempts to used advanced functions to fix the cross correlation, 
#			this was slower and didn't work
#	5	-	Most of the individual functions broken out into a intf_tools library 
#			the global settings object introduced
#	6	-	Major revision to the file format, now including a header.  
#			Added capability for non perpendicular baselines, and more than 3 antennas
#	6.1	-	Moved the intf_process argument parser out of intf_tools (why was it there)


version = 6.1

import sys,os,time,argparse
import struct

#these libraries may not exist, but are needed
try:
	from numpy import *
except:
	print 'ERROR: numpy does not appear to be installed!'
	sys.exit(2)



#this sets the orders for the antennas, it should probably be in the
#Settings object
AntOrders = [ 	[1, 0, 2],
				[3, 2, 0],
				[0, 3, 1],
				[2, 1, 3]  ]

#########
#Helper Functions
import intf_tools6 as st

def pkpk(arr):
	return arr.max() - arr.min()

def rms(arr):
	return (mean(arr**2))**.5

##########
#Main Functions

def process_window(data,lastEvent,ants,Settings):
	"""
	data      - time waveforms to be correlated
	lastEvent - Needed to determine dups, can be None
	ants      - array of the index of the antennas to use, ordered
	Settings  - Settings object with global settings
	"""
	N  = Settings.numSamples
	I  = Settings.numIterate
	P  = Settings.numInterp
	bw = Settings.bandwidth
	
	basel = st.antBaseline( ants )
	theta = st.antTheta( ants )
	dtheta = theta[0]-theta[1]
	
	isNoise = False
	isDup   = False	
	
	Wt = sin( arange( N ) *pi /N )*2./3+1./3
	
	###
	#Step 1: Initial prepping of the data
	
	#calculate the time weighting, it's a modified hamming window
	#modified by being raised higher
	Wt_raise = 0.16	#for hamming, Wt_raise = 0.08
	Wt = (1-cos( 2*pi*arange(N)/N ))*(1-Wt_raise)*.5+Wt_raise
	
	#weight the data towards the center of the window
	d0 = ( data[ants[0],:] - mean(data[ants[0],:]) )*Wt
	d1 = ( data[ants[1],:] - mean(data[ants[1],:]) )*Wt
	d2 = ( data[ants[2],:] - mean(data[ants[2],:]) )*Wt	

	#calculate the ffts
	fft0 = fft.fft(d0,2*N)
	fft1 = fft.fft(d1,2*N)
	fft2 = fft.fft(d2,2*N)
	
	#this is used for duplicate testing, 
	#d0 is the center of the array, so timing is based on this channel
	iMax = argmax( abs(d0) ) + Settings.startSample	

	#check for a duplicate event
	if lastEvent != None:		
		if 	abs(lastEvent.iMax - iMax) < Settings.minSep:
			#it's a duplicate
			isDup = True	
		elif ants != AntOrders[0]:
			#we can't make a new event yet!
			isDup   = True
			isNoise = True
	else:
		if ants != AntOrders[0]:
			return None

	###
	# Step 2: Calculate the windowing function
	freq = abs( 180*fft.fftfreq( 2*N ) )
	W = zeros(2*N)+1
	#there's a very very strong transmitter belower 20 MHz, it
	#blows through the filter pretty well and needs to be hit pretty 
	#hard to remove it
	W[ freq<bw[0] ] = 0	#The stuff outside the hardware filters
	W[ freq>bw[1] ] = 0

	#calculate some metrics now while we have a simple windowing function
	#calculate the autospectra of the signals
	d0 = real( fft.ifft( fft0*W ) )
	d1 = real( fft.ifft( fft1*W ) )
	d2 = real( fft.ifft( fft2*W ) )
	#0 offset is index 0
	eamp0 = sqrt( min( pkpk(d0),pkpk(d1),pkpk(d2) ) )
	eamp1 = sqrt( min( rms(d0) ,rms(d1) ,rms(d2)  ) )
	
	#calculate the color
	cSig = 6.25	#standard deviation of the color filters
	cRed = exp(-.5*( (freq-42.5)/cSig )**2)
	cGrn = exp(-.5*( (freq-55.0)/cSig )**2)
	cBlu = exp(-.5*( (freq-67.5)/cSig )**2)
	#integrate
	cRed = sqrt(sum(cRed*abs(fft0*W)**2))
	cGrn = sqrt(sum(cGrn*abs(fft0*W)**2))
	cBlu = sqrt(sum(cBlu*abs(fft0*W)**2))
	#normalize
	cMax = max(cRed,cGrn,cBlu)
	cRed = int( cRed/cMax*255 )
	cGrn = int( cGrn/cMax*255 )
	cBlu = int( cBlu/cMax*255 )

	
	W0 = W.copy()	#used for a small number of calculations which require a simple windowing function
	
	# use maximum likelihood style filtering?
	if Settings.windowing == None or Settings.windowing == 0:
		#default case
		W01 = W02 = W12 = W
	#the old ML window
	if Settings.windowing == 1:
		#using the WII filter (I hope)
		#Estimate the signal spectra
		key = tuple( sorted( (ants[0],ants[1]) ))
		Y01 = st.symm_smooth( abs(fft0*fft1.conj()), 2 ) - \
				abs( st.Spectra[ key ] )
		Y01[Y01<0] = 0
		key = tuple( sorted( (ants[0],ants[2]) ))
		Y02 = st.symm_smooth( abs(fft0*fft2.conj()), 2 ) - \
				abs( st.Spectra[ key ] )
		Y02[Y02<0] = 0
		key = tuple( sorted( (ants[1],ants[2]) ))
		Y12 = st.symm_smooth( abs(fft1*fft2.conj()), 2 ) - \
				abs( st.Spectra[ key ] )
		Y12[Y12<0] = 0

		#load the noise auto spectra	
		AN0 = abs(st.Spectra[ants[0]])
		AN1 = abs(st.Spectra[ants[1]])
		AN2 = abs(st.Spectra[ants[2]])

		W01 = Y01/(AN0*AN1 + Y01*(AN0+AN1) + Y01**2)
		W01[ (freq<bw[0])|(freq>bw[1]) ] = 0
		#normalize
		if not max(W01) == 0:
			W01 /= max(W01)

		W02 = Y02/(AN0*AN2 + Y02*(AN0+AN2) + Y02**2)
		W02[ (freq<bw[0])|(freq>bw[1]) ] = 0
		#normalize
		if not max(W02) == 0:
			W02 /= max(W02)

		W12 = Y12/(AN1*AN2 + Y12*(AN1+AN2) + Y12**2)
		W12[ (freq<bw[0])|(freq>bw[1]) ] = 0
		#normalize
		if not max(W12) == 0:
			W12 /= max(W12)
		
	if Settings.windowing == 2:	#PHAT
		nSmooth = 5
		#Gs_00 = st.symm_smooth( abs(fft0)**2, nSmooth )
		#Gs_00[0] += max(Gs_00)
		#Gs_11 = st.symm_smooth( abs(fft1)**2, nSmooth )
		#Gs_11[0] += max(Gs_11)
		#Gs_22 = st.symm_smooth( abs(fft2)**2, nSmooth )
		#Gs_22[0] += max(Gs_22)		
		Gs_01 = st.symm_smooth( fft0*fft1.conj(), nSmooth )
		Gs_01[0] += max(Gs_01)
		Gs_02 = st.symm_smooth( fft0*fft2.conj(), nSmooth )
		Gs_02[0] += max(Gs_02)
		Gs_12 = st.symm_smooth( fft1*fft2.conj(), nSmooth )
		Gs_12[0] += max(Gs_12)
		
		W01 = 1./abs( Gs_01 )
		W02 = 1./abs( Gs_02 )
		W12 = 1./abs( Gs_12 )
	
	if Settings.windowing == 3: #SCOT
		nSmooth = 5
		Gs_00 = st.symm_smooth( abs(fft0)**2, nSmooth )
		Gs_00[0] += max(Gs_00)
		Gs_11 = st.symm_smooth( abs(fft1)**2, nSmooth )
		Gs_11[0] += max(Gs_11)
		Gs_22 = st.symm_smooth( abs(fft2)**2, nSmooth )
		Gs_22[0] += max(Gs_22)		
		#Gs_01 = st.symm_smooth( fft0*fft1.conj(), nSmooth )
		#Gs_01[0] += max(Gs_01)
		#Gs_02 = st.symm_smooth( fft0*fft2.conj(), nSmooth )
		#Gs_02[0] += max(Gs_02)
		#Gs_12 = st.symm_smooth( fft1*fft2.conj(), nSmooth )
		#Gs_12[0] += max(Gs_12)
		
		W01 = 1./sqrt( Gs_00 * Gs_11 )
		W02 = 1./sqrt( Gs_00 * Gs_22 )
		W12 = 1./sqrt( Gs_11 * Gs_22 )
	
	if Settings.windowing == 4: #ML
	        nSmooth = 5
		Gs_00 = st.symm_smooth( abs(fft0)**2, nSmooth )
		Gs_00[0] += max(Gs_00)
		Gs_11 = st.symm_smooth( abs(fft1)**2, nSmooth )
		Gs_11[0] += max(Gs_11)
		Gs_22 = st.symm_smooth( abs(fft2)**2, nSmooth )
		Gs_22[0] += max(Gs_22)		
		Gs_01 = st.symm_smooth( fft0*fft1.conj(), nSmooth )
		Gs_01[0] += max(Gs_01)
		Gs_02 = st.symm_smooth( fft0*fft2.conj(), nSmooth )
		Gs_02[0] += max(Gs_02)
		Gs_12 = st.symm_smooth( fft1*fft2.conj(), nSmooth )
		Gs_12[0] += max(Gs_12)		
		
		Gam_01 = abs( Gs_01/( Gs_00*Gs_11 )**.5 )
		Gam_01[ Gam_01 > 1] = 1
		Gam_02 = abs( Gs_02/( Gs_00*Gs_22 )**.5 )
		Gam_02[ Gam_01 > 1] = 1
		Gam_12 = abs( Gs_12/( Gs_11*Gs_22 )**.5 )
		Gam_12[ Gam_01 > 1] = 1
		
		W01 = abs(Gam_01)**2/( abs(Gs_01)*(1+1e-6-abs(Gam_01)**2))
		W02 = abs(Gam_02)**2/( abs(Gs_02)*(1+1e-6-abs(Gam_02)**2))
		W12 = abs(Gam_12)**2/( abs(Gs_12)*(1+1e-6-abs(Gam_12)**2))
		
		

	
	###
	#Step 3: Calculate the Correlation
        
        #We need to recalculate the autospectra to normalize things
	a01 =  real( fft.ifft( fft0*W01*fft0.conj() ) ) *\
	       real( fft.ifft( fft1*W01*fft1.conj() ) )
	a02 =  real( fft.ifft( fft0*W02*fft0.conj() ) ) *\
	       real( fft.ifft( fft2*W02*fft2.conj() ) )
	
	#calculate the cross spectrum
	x01 = W01 * fft0 * (fft1.conj())
	x02 = W02 * fft0 * (fft2.conj())
	x12 = W12 * fft1 * (fft2.conj())
	#zeropad the cross spectrum so the xcorr is interpolated
	x01 = st.fpad(x01,P)
	x02 = st.fpad(x02,P)
	x12 = st.fpad(x12,P)
	#invert and normalize
	x01 = fft.fftshift( fft.ifft(x01) )[1:]
	x02 = fft.fftshift( fft.ifft(x02) )[1:]
	x12 = fft.fftshift( fft.ifft(x12) )[1:]
	#normalization is window function dependant
	x01 /= sqrt(a01[0])
	x02 /= sqrt(a02[0])
	#drop the imaginary part
	x01 = real( x01 )
	x02 = real( x02 )
	x12 = real( x12 )	
	
	###
	# Step 4: Calculate location
	
	#s calculated in us
	s01 = (st.xmax(x01)/P -N)/Settings.sampleRate
	s02 = (st.xmax(x02)/P -N)/Settings.sampleRate
	s12 = (st.xmax(x12)/P -N)/Settings.sampleRate
	#offset by the cable delays
	d01 = (Settings.antDels[ants[0]] - Settings.antDels[ants[1]])/1000.
	d02 = (Settings.antDels[ants[0]] - Settings.antDels[ants[2]])/1000.
	d12 = (Settings.antDels[ants[1]] - Settings.antDels[ants[2]])/1000.
	
	s01 = s01 - d01
	s02 = s02 - d02
	s12 = s12 - d12
	
	#more metrics
	expeak = min( [x01.max(),x02.max()] )
	eclos  = (s01+s12-s02)*1000
	efwhm  = 0	#not used anymore
	
	#the angles of incidence
	#print basel
	cosa = s01/1e6*Settings.c/ basel[0]
	cosb = s02/1e6*Settings.c/ basel[1]
	#print st.cosab2AzEl(cosa, cosb, ants)
	#convert cosa and cosb into perpendular coords
	cosa, cosb = st.cosab2cosab( cosa, cosb, ants )
	#print cosa**2+cosb**2
	#print st.cosab2AzEl_perp(cosa, cosb)

	#is the solution valid?
	if cosa**2 + cosb**2 > 1.1:
		#~ #do the noise stuff
		if Settings.windowing == None or Settings.noise:
			st.store_spectra( [fft0,fft1,fft2], ants)

		return None

	if isDup:
		if not isNoise:
			lastEvent.append( eamp0,eamp1,expeak,eclos,cosa,cosb )

		#~ print 'appending old event',Settings.startSample
		return lastEvent

		
	event = st.Event(Settings.startSample,eamp0,eamp1,expeak,eclos,cosa,cosb,ants,iMax, (cRed,cGrn,cBlu))
	#~ print 'returning new event',Settings.startSample
	return event		



def parse_args():

	parser = argparse.ArgumentParser(description="Process DITF data",fromfile_prefix_chars='@')
	parser.add_argument('--version', action='version', version='%(prog)s'+' %f'%version)
	
	#the implicit arguments
	parser.add_argument('input_file',
		help='Path to binary input file')
	parser.add_argument('output_file', nargs='?', default=None,
		help='Path to ascii output file')
	
	#the primary arguments
	parser.add_argument('-S','--numSamples',default=256,type=int,
		help='The number of samples in a window')
	parser.add_argument('-I','--numIterate',default=64,type=int,
		help='The number of samples to incriment between windows, if this is smaller than numSamples the windows will overlap')
	parser.add_argument('-P','--numInterp',default=2,type=int,
		help='The factor to upsample by, 16 increases the number of samples by a factor of 16')
	parser.add_argument('-s','--start',default=0,type=float,
		help='Start position in ms from the begininng of the file')
	parser.add_argument('-p','--stop',default=-1,type=float,
		help='Stop position in ms from the beginning of the file')
	parser.add_argument('-W','--windowing',default=None,type=int,
		help='Enable speciial windowing to enhance the correlation')
	parser.add_argument('-A','--ants',default=3,type=int,
		help='The number of antennas to use for the solutions')
	parser.add_argument('-N','--noise_file',default=None,type=str,
		help='Path to the xns file containing the noise specta')

	#Debugging and verbose options
	parser.add_argument('-v','--verbose',action='count',
		help='Print output to std_out as well as to file')
	parser.add_argument('--noise',action='store_true',default=False,
		help='This assumes all solutions are noise, useful for getting the noise spectra only')

	#parse the command line
	arguments = parser.parse_args(sys.argv[1:])
	st.Settings.arguments = arguments

	#input output files
	st.Settings.inFileS = arguments.input_file
	if not os.path.exists(st.Settings.inFileS):
		print 'input file does not exist: %s'%st.Settings.inFileS
		sys.exit(1)
	if arguments.output_file != None:
		st.Settings.outFileS = arguments.output_file
	else:
		st.Settings.outFileS = os.path.splitext(st.Settings.inFileS)[0]+'.dat.gz'

	if arguments.noise_file != None:
		st.Settings.noiseFileS = arguments.noise_file
	else:
		st.Settings.noiseFileS = os.path.splitext(st.Settings.outFileS)[0]+'.xns'

	#set the processing characteristics
	st.Settings.numSamples = arguments.numSamples
	st.Settings.numIterate = arguments.numIterate
	st.Settings.numInterp  = arguments.numInterp
	st.Settings.windowing = arguments.windowing
	st.Settings.minSep    = st.Settings.numSamples/2
	st.Settings.numAnts   = arguments.ants
	st.Settings.verbose   = arguments.verbose
	st.Settings.noise     = arguments.noise


	#read the input file header and store that
	st.Settings.header = st.RawHeader(st.Settings.inFileS)
	
	#from the header, we want to have the sample rate and pretrigger length on hand
	st.Settings.sampleRate  	  = st.Settings.header.sampleRate/1000000	#in samples/us
	st.Settings.preTriggerSamples = st.Settings.header.preTriggerSamples
	
	#deal with the start and stop samples, surprising complicated, this
	if arguments.start == 0:
		#that's the beginning of the file, easy
		st.Settings.startSample = 0
	else:
		st.Settings.startSample = \
			arguments.start*st.Settings.sampleRate*1000+st.Settings.preTriggerSamples
	if arguments.stop == -1:	#process the whole file
		st.Settings.stopSample = -1
	else:
		st.Settings.stopSample = \
			arguments.stop*st.Settings.sampleRate*1000+st.Settings.preTriggerSamples


	#Are we saving or loading the noise spectra?
	if st.Settings.windowing != None and not st.Settings.noise and st.Settings.windowing != 0:
		st.Settings.storeSpec = False
	else:
		st.Settings.storeSpec = True
	
	#load the noise spectra?
	if not st.Settings.storeSpec and st.Settings.windowing != 0:
		st.read_noise_spectra()
	else:
		st.init_noise_spectra()
		
	
def main():
	#get the arguments
	#build the settings
	parse_args()
	
	if st.Settings.verbose:
		print st.Settings.header.report()
		print st.Settings.report()
	
	st.Settings.file_initialize()
	
	global AntOrders
	if st.Settings.numAnts == 3:
		#there's 1 permutation
		AntOrders = AntOrders[:1]
	elif st.Settings.numAnts == 4:
		#there's 4 permutations
		AntOrders = AntOrders[:4]
	
	lastEvent = None
	saturation_counter = 0
	
	#remember to remove this later!!!
	#st.Settings.startSample = 702389312
	
	while st.Settings.startSample <= st.Settings.stopSample or st.Settings.stopSample == -1:
		#for noise processing, this is the progress bar
		if st.Settings.noise:
			print st.Settings.stopSample-st.Settings.startSample

		st.Settings.startSample += st.Settings.numIterate

		###
		#Step 1: Read data
		data = st.read_raw_file_data()
		if data == None:
			#we're out of data
			break
		
		#check for saturation on all channels
		#~ if data.max()==2**16 or data.min()==0:
			#~ print '**** Saturation Detected ****'
			#~ st.Settings.SatCounter += 1
			#~ #move forward farther than normal
			#~ #st.Settings.startSample += st.Settings.numSamples/2
			#~ continue
		
		#process the window
		for order in AntOrders:
			event = process_window(data,lastEvent, order, st.Settings)
			#event handling
			if event == None:
				continue
			if lastEvent != None:
				if not lastEvent.iMax == event.iMax:
					#write the event to the output file
					lastEvent.write(st.Settings.outFileS)
					if st.Settings.verbose >= 1:
						lastEvent.printout()
						#~ #print len(lastEvent.iClos)
						#~ print lastEvent.ants, len(lastEvent.iClos)
						#~ for i in range( len(lastEvent.iClos) ):
							#~ print lastEvent.cosa[i][0], lastEvent.cosb[i][0], lastEvent.iClos[i]
						#print lastEvent.cosa
						#print lastEvent.cosb
						#sys.exit()
			lastEvent = event
			

		#sys.exit()
	
	#the main loop is done, finalize the file
	st.Settings.file_finalize()
	
	#check if we've processed all requested data
	#store the spectrum info
	if st.Settings.windowing == None or st.Settings.noise:
		st.write_noise_spectra()
		
		
	
if __name__ == '__main__':
	main()

		
	
	
