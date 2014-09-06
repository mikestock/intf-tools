#!/usr/bin/python
#######
# intf_process_multi.py <FILE_IN> [FILE_OUT]
#
#	Reads FILE_IN and calculates azimuth/elevation for the data
#	Input is AlazarTech 9642 4 channel digitized data, should be 
#	version 12+
#	If FILE_OUT is supplied, writes locations and metrics.  
#
#	Call with --help for information about arguments
#
# Written by Michael Stock, 2014
# versions:
#	1	-	Initial code, sort of worked
#	2.0	-	Transitioned to a fully multiprocess enabled code
#	2.1	-	Added a read buffer, which helps greatly when used with 
#			many threads

version = 2.1

import sys,os,time,argparse, struct
import multiprocessing as mp

#these libraries may not exist, but are needed
try:
	from numpy import *
except:
	print 'ERROR: numpy does not appear to be installed!'
	sys.exit(2)

import intf_tools as st
it = st	#match the new and old way to referencing this
import cintf_tools as cit

#this sets the orders for the antennas, it should probably be in the
#Settings object
#2013 data

#########
# Helper Functions
def pkpk(arr):
	return arr.max() - arr.min()

def rms(arr):
	return (mean(arr**2))**.5

def find_peak( X, x0):
	thresh = max(X)/10
	###
	# written for easy implementation in Cython
	
	#Left Peak
	iML = x0
	mML = X[x0]
	#Right Peak
	iMR = x0
	mMR = X[x0]

	RLoop2 = True

	###
	#Left Loop part 1
	i = x0-1
	#print 'L loop 1'
	while X[i] > thresh:
		#print i, 'guess pos'
		if X[i] > mML:
			iML = i
			mML = X[i]
		i -= 1
		if i < 0:
			break
	#print 'L max0', iML, x0
	#Left Loop part 2
	if iML == x0 and i > 0:
		if X[x0] > thresh:
			RLoop2 = False
		mML = 0
		#print 'L loop 2'
		while X[i] < thresh:
			#print i, 'guess neg'
			i -= 1
			if i < 0:
				break
		if i >= 0:
			#loop down till X goes negaitve
			while X[i] > 0:
				#print i, 'guess pos'
				if X[i] > mML:
					iML = i
					mML = X[i]
				i -= 1
				if i < 0:
					break
	#print 'L Max',iML
		
	###
	#Right Loop part 1
	i = x0+1
	#print 'R loop 1'
	while X[i] < thresh:
		#print i, 'guess pos'
		if X[i] > mMR:
			iMR = i
			mMR = X[i]
		i += 1		
		if i >= X.shape[0]:
			break

	if i < X.shape[0]:
		#then we still have array to search
		while X[i] > thresh:
			#print i, 'guess pos'
			if X[i] > mMR:
				iMR = i
				mMR = X[i]
			i += 1
			if i >= X.shape[0]:
				break
	else:
		#well fuck, we're still in the negative part of the 
		#xc, and didn't find what we were looking for.  
		#this should result is the right loop not producing 
		#a solution
		pass
	#print 'R max', iMR, x0
	#Right Loop part 2
	#We don't want to do this if we started on a hill
	if iMR == x0 and RLoop2 and i < X.shape[0]:
		mMR = 0
		#print 'R loop 2'
		while X[i] < thresh:
			i += 1
			if i >= X.shape[0]:
				break
		#loop down till X goes negaitve
		if i < X.shape[0]:
			while X[i] > thresh:
				#print i, 'guess pos'
				if X[i] > mMR:
					#print i, 'new max'
					iMR = i
					mMR = X[i]
				i += 1
				if i >= X.shape[0]:
					break
	#print 'R Max',iMR

	#determine which peak, right or left, is closer to the initial guess
	if x0-iML < iMR - x0:
		iMax = iML
	elif iMR - x0 < x0-iML:
		iMax = iMR
	else:	#they're equal, damnit
		if mML > mMR:
			iMax = iML
		else:
			iMax = iMR
	
	#print 'iMax',iMax
	
	#now that all of that is done, we make a parabolic fit to the 
	#points around iMax
	#this used to be done wit built in functions, but they're apparently 
	#slow
	
	###
	# step 1, get the slops (reusing varables here)
	mML = X[iMax] - X[iMax-1]	#value taken at -1/2 (in samples)
	mMR = X[iMax+1] - X[iMax]	#value taken at +1/2 (in samples)
	#solving now for y = mx + b
	b = (mML + mMR)/2	#the +/- .5 will always cancel the m's
	m = mMR - mML	#the b's here alway cancel, and the .5's add to 1
	#and, this is the result!
	#return iML, iMR, iMax-b/m, thresh
	return iMax-b/m

#########
# Multiprocessing Classes
class WindowProcessor( mp.Process):
	
	def __init__( self, inQ, outQ, noiseQ, M, Settings ):
		
		#initialize mp stuff
		mp.Process.__init__(self)
		
		#initialize intf stuff
		self.inQ  = inQ		#where we get raw data
		self.outQ = outQ	#where the processed data goes
		self.noiseQ=noiseQ	#for noise spectra computation
		self.M    = M		
		self.Settings = Settings
		self.running  = True
	
	def run(self):
		######
		# initialize
		Settings = self.Settings
		M  = self.M
		S  = Settings.numSamples
		I  = Settings.numIterate
		P  = Settings.numInterp
		bw = Settings.bandwidth
		pix= Settings.pix
		windowing = Settings.windowing
		nSmooth = 5
		fs = 180.

		#calculate the time weighting, it's a modified hamming window
		#modified by being raised higher
		Wt_raise = 0.16	#for hamming, Wt_raise = 0.08
		Wt = (1-cos( 2*pi*arange(S)/S ))*(1-Wt_raise)*.5+Wt_raise
		
		#this is used for the centroid calculations
		ica, icb = mgrid[ slice(pix), slice(pix) ]
		#this is the analog interferometer resolution
		ares = 300./bw[1]/Settings.arrDiam

		
		while self.running:
			######
			# Get data from the queue
			sSample, data = self.inQ.get()

			event = None

			######
			# Generate the windowing function
			#generate the simplistic windowing function to bandwidth limit 
			#the cross correlations
			W = zeros(2*S)+1
			freq = abs( 180*fft.fftfreq( 2*S )  )
			W[ abs(freq) < bw[0] ] = 0
			W[ abs(freq) > bw[1] ] = 0


			######
			# Compute the FFT's
			ffts  = empty( (M, 2*S), dtype='complex' )
			erms  = zeros( M )
			epkpk = zeros( M )
			for i in range(M):
				d = (data[i,:] - mean(data[i,:]))*Wt
				ffts[i,:] = fft.fft( d, 2*S )
				#we track the peak of the 0 channel for overlapping windows
				if i == 0:
					iMax = argmax( abs(d) ) + sSample
				#track the amplitude
				d = real(fft.ifft( ffts[i,:]*W )[:S])
				epkpk[i] =  st.pkpk(d)
				erms[i] =  st.rms(d) 

			###
			# do a test on the rms values, 
			# they should be all about the same
			if std(erms) > mean(erms):
				#don't do noise stuff, the window is badly contaminated
				self.outQ.put( [sSample, iMax, event] )
				continue				

			# Compute the crossspectra and Signal Spectra
			G   = zeros( [M*(M-1)/2, S*2], dtype='complex' )		#the cross spectra
			Wij = zeros( [M*(M-1)/2, S*2] )		#the weighting function
			Y   = zeros( [S*2] )	#the signal spectra
			k = 0
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					G [ k,:] = ffts[i,:] * (ffts[j,:].conj())
					
					if windowing in [0,2,3,None]:
						#that's none, PHAT, or SCOT
						Y[:] = 1.
					elif i in [0,1,2] and j in [0,1,2]:
						#estimate the signal spectra
						#omit antenna 3 (it's noisy)
						Y += abs(G[k])/3
					
					k += 1
			
			######
			# Correlate
			A    = zeros( [M*(M-1)/2, 2] )	    #array, baseline angle
			xcs  = zeros( [M*(M-1)/2, S*2*P] )	#store the xcross correlations in an array
			dls  = zeros( M*(M-1)/2 )	        #store the delays in an array
			bls  = zeros( M*(M-1)/2 )	        #sotre the baselines in an array
			k = 0		#this is used to track the index of the arrays
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					th = it.antTheta([i,j])[0]
					bl = it.antBaseline([i,j])[0]
					
					###
					# Calculate the windowing function

					if windowing == None or windowing == 0:
						#no windowing
						Wij = W
					elif windowing == 1:
						#WII windowing

						#this is the old code
						#~ Y12 = st.symm_smooth( abs(fft1*fft2.conj()), 2 ) - \
							#~ abs( st.NoiseSpectra[ key ] )
						#~ Y12[Y12<0] = 0

						# this one is complicated, and makes use of the
						# signal spectra Y
						ANi = abs(st.NoiseSpectra[i])
						ANj = abs(st.NoiseSpectra[j])
						ANij= abs(st.NoiseSpectra[(i,j)])
						
						# Get the signal spectra for this window 
						# and subtract off known correlated noise
						y      = Y - ANij
						#~ y = it.symm_smooth( G[k], 2 ) - ANij
						y[y<0] = 0
						
						
						# this is the windowing function
						Wij = W*y/( ANi*ANj + y*(ANi+ANj) + y*y )
						# Normalize to 1, there's not a good reason 
						# to do this really.  
						if not max(Wij) == 0:	#don't divide by 0
							Wij /= max(Wij)
					elif windowing == 2:
						#PHAT
						# simply divide by the amplitude
						w = abs(G[k])
						w[w<1e-10] = 1e-10 #don't divide by 0
						Wij = W/w
					elif windowing == 3:
						#SCOT
						# uses a smooth auto correlation
						w  = it.symm_smooth( abs(ffts[i])**2, nSmooth )
						w *= it.symm_smooth( abs(ffts[j])**2, nSmooth )
						w[w<1e-10] = 1e-10 #don't divide by 0
						Wij = W/sqrt(w)
					else:
						#I don't know what this is
						print 'Unknown Windowing Option'
						sys.exit(windowing)
					
					###
					# interpolate
					g = it.fpad(G[k]*Wij, P)
					
					###
					# the cross correlation
					X =  real(fft.fftshift( fft.ifft(g) ))
					
					###
					# Normalization
					# this won't affect the peak, but does affect correlation amplitude
					# note, expensive way to calculate this
					aci = real( fft.ifft( ffts[i,:]*Wij*ffts[i,:].conj() ) )[0]
					acj = real( fft.ifft( ffts[j,:]*Wij*ffts[j,:].conj() ) )[0]
					if aci > 0 and acj > 0:
						X /= sqrt(aci*acj)
					#there's a factor of 1.5 difference between the correlation values 
					#in the 2 baseline algorithm, and here.  
					#I have not been able to track it down.  For now, so they noise 
					#reduction stays about the same for both, just correct for it.
					X  = real(X)
					###
					# Store the arrays
					xcs[k,:] = X
					bls[k]   = bl
					dls[k]   = (Settings.antDels[i] - Settings.antDels[j])/1000
					A[k,0]   = cos(th)
					A[k,1]   = sin(th)
					
					k += 1
			
			######
			# Initial Guess with Imaging
			# the pseudo image routine requires arrays as input
			###
			# compute the pseudo image
			bMask   = ones(xcs.shape[0], dtype='bool')
			
			
			bbox = array( [[-1.5,1.5],[-1.5,1.5]] )
			h= cit.proj_image(xcs[bMask,:], 
							  bls[bMask], 
							  dls[bMask], 
							  A[bMask,:], N=pix, fs=180*P, bbox=bbox)
			hmax = unravel_index(h.argmax(), h.shape)					
			
			###
			# The initial guess
			cosa = (bbox[0,1]-bbox[0,0])*(hmax[0]+.5)/pix+bbox[0,0]
			cosb = (bbox[1,1]-bbox[1,0])*(hmax[1]+.5)/pix+bbox[1,0]

						
			###
			# Find baselines far from the source
			###
			# First, compute the time delays
			mMiddle = (xcs.shape[1]+1)/2
			expeak  = zeros( A.shape[0] )
			taus = {}
			
			k = 0
			iR = xcs.shape[0]
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue

					#check for a record with no signal!
					if xcs[k].max() == 0:
						expeak[k] = 0
						taus[(i,j)] = 0,k
						iR -= 1
						k += 1
						continue

					###
					# Calculate time delay
					tau = (A[k,0]*cosa+A[k,1]*cosb)*bls[k]/300. - dls[k]
					#convert to samples
					m0   = mMiddle - round(tau*180*P)

					###
					# parabolic fit to nearest peak
					# note, this operation is slow
					#m   = cit.find_peak(xcs[i], m)
					m   = find_peak(xcs[k], m0)

					###
					# should we keep this baseline?
					#print ares*bls[k]/300.*fs*P/3
					#print m0-m
					if 	abs(m0-m) > ares*bls[k]/300.*fs*P/2 or \
						bMask[k] == False:
						#this solution is far off base
						bMask[k] = False
						expeak[k]   = 0
						k+=1
						taus[(i,j)] = 0,-1
						iR -= 1
						continue

					expeak[k] = xcs[k, round(m) ]
					#Calculate the time delay for closure
					sij  = (m/P-S)/180.
					sij -= dls[k]
					sij  = -sij
					taus[(i,j)] = sij*1000, k

					k += 1
			if iR < 2:
				#require at least 2 of the local peaks to occur reasonably 
				#close the the global peak
				self.outQ.put( [sSample, iMax, event] )
				#~ print 'removing', 453
				continue

			######
			# Second Guess with Imaging
			# the pseudo image routine requires arrays as input
			###
			# compute the pseudo image
			bbox = array( [[cosa-ares, cosa+ares],[cosb-ares,cosb+ares]] )
			#~ h= cit.proj_image(xcs, 
							  #~ bls, 
							  #~ dls, 
							  #~ A, N=pix, fs=180*P, bbox=bbox)
			#print bMask
			h= cit.proj_image(xcs[bMask,:], 
							  bls[bMask], 
							  dls[bMask], 
							  A[bMask,:], N=pix, fs=180*P, bbox=bbox)

			###
			# Find the maxima
			hmax = unravel_index(h.argmax(), h.shape)
			
			##
			#remove all values over 1 ares away
			#this is the size of a pixel
			dpix = (bbox[0,1]-bbox[0,0])/pix
			r = ( (ica-ica[hmax[0],hmax[1]])**2 + (icb-icb[hmax[0],hmax[1]])**2 )
			h[ r>(ares/dpix)**2 ] = 0
			##
			#remove all values lower than the half max
			h[ h<.5*h.max() ] = 0

			#compute the centroid (in pixels)
			cosa = sum( h*ica )/ sum(h)
			cosb = sum( h*icb )/ sum(h)
			#convert to sky coords
			cosa = (bbox[0,1]-bbox[0,0])*(cosa+.5)/pix+bbox[0,0]
			cosb = (bbox[1,1]-bbox[1,0])*(cosb+.5)/pix+bbox[1,0]

			###
			# is the solution valid? (almost imposible not to be at this point)
			if cosa**2 + cosb**2 > 1.1:
				#~ print 'removing', 496
				#do the noise stuff
				if Settings.windowing == None or Settings.noise:
					self.noiseQ.put( ffts )

				self.outQ.put( [sSample, iMax, event] )
				continue
					

			###
			# reduce the metrics
			# the :3 is here for compatibility reasons
			erms = mean(erms[:3] )
			epkpk = mean(epkpk[:3] )

			#####
			# Closure
			# this is a pain to get now


			eclos = []
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					for k in range(j+1,M):
						if j == k:
							continue
						#check that all combinations are in eMask 
						#(included in the solution)
						if  taus[i,j][1] >= 0 and \
							taus[j,k][1] >= 0 and \
							taus[i,k][1] >= 0:
							eclos.append( taus[i,j][0]+taus[j,k][0]-taus[i,k][0] )
			
			###
			# This filters out badly damaged solutions
			if len(eclos) == 0 :#or max(eclos) > 4:
				#~ print 'removing', 534
				#~ print 'closure error'
				#do the noise stuff
				if Settings.windowing == None or Settings.noise:
					self.noiseQ.put( ffts )

				self.outQ.put( [sSample, iMax, event] )
				continue
				
			eclos = min(eclos)
			
			###
			# The correlation peak has been problematic
			# The 0,1,2 antennas correlate better than the 3 antenna
			# so, omit correlations with this antenna
			mask = zeros( M*(M-1)/2, dtype='bool')
			k = 0
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					if i in [0,1,2] and j in [0,1,2]:
						mask[k] = True
					k += 1
			expeak= mean(expeak[mask])

			event = [erms,epkpk,expeak,eclos,cosa,cosb,iR,(iR,0,0)]
			self.outQ.put( [sSample, iMax, event] )

	def run_bakc(self):
		######
		# initialize
		Settings = self.Settings
		M  = self.M
		S  = Settings.numSamples
		I  = Settings.numIterate
		P  = Settings.numInterp
		bw = Settings.bandwidth
		windowing = Settings.windowing
		nSmooth = 5

		#calculate the time weighting, it's a modified hamming window
		#modified by being raised higher
		Wt_raise = 0.16	#for hamming, Wt_raise = 0.08
		Wt = (1-cos( 2*pi*arange(S)/S ))*(1-Wt_raise)*.5+Wt_raise

		while self.running:
			######
			# Get data from the queue
			sSample, data = self.inQ.get()

			event = None

			######
			# Generate the windowing function
			#generate the simplistic windowing function to bandwidth limit 
			#the cross correlations
			W = zeros(2*S)+1
			freq = abs( 180*fft.fftfreq( 2*S )  )
			W[ abs(freq) < bw[0] ] = 0
			W[ abs(freq) > bw[1] ] = 0


			######
			# Compute the FFT's
			ffts  = empty( (M, 2*S), dtype='complex' )
			erms  = zeros( M )
			epkpk = zeros( M )
			for i in range(M):
				d = (data[i,:] - mean(data[i,:]))*Wt
				ffts[i,:] = fft.fft( d, 2*S )
				#we track the peak of the 0 channel for overlapping windows
				if i == 0:
					iMax = argmax( abs(d) ) + sSample
				#track the amplitude
				d = real(fft.ifft( ffts[i,:]*W )[:S])
				epkpk[i] =  st.pkpk(d)
				erms[i] =  st.rms(d) 

			# Compute the crossspectra and Signal Spectra
			G   = zeros( [M*(M-1)/2, S*2], dtype='complex' )		#the cross spectra
			Wij = zeros( [M*(M-1)/2, S*2] )		#the weighting function
			Y   = zeros( [S*2] )	#the signal spectra
			k = 0
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					G [ k,:] = ffts[i,:] * (ffts[j,:].conj())
					
					if windowing in [0,2,3,None]:
						#that's none, PHAT, or SCOT
						Y[:] = 1.
					elif i in [0,1,2] and j in [0,1,2]:
						#estimate the signal spectra
						#omit antenna 3 (it's noisy)
						Y += abs(G[k])/3
					
					k += 1
			
			######
			# Correlate
			A    = zeros( [M*(M-1)/2, 2] )	    #array, baseline angle
			xcs  = zeros( [M*(M-1)/2, S*2*P] )	#store the xcross correlations in an array
			dls  = zeros( M*(M-1)/2 )	        #store the delays in an array
			bls  = zeros( M*(M-1)/2 )	        #sotre the baselines in an array
			k = 0		#this is used to track the index of the arrays
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					th = it.antTheta([i,j])[0]
					bl = it.antBaseline([i,j])[0]
					
					###
					# Calculate the windowing function

					if windowing == None or windowing == 0:
						#no windowing
						Wij = W
					elif windowing == 1:
						#WII windowing

						#this is the old code
						#~ Y12 = st.symm_smooth( abs(fft1*fft2.conj()), 2 ) - \
							#~ abs( st.NoiseSpectra[ key ] )
						#~ Y12[Y12<0] = 0

						# this one is complicated, and makes use of the
						# signal spectra Y
						ANi = abs(st.NoiseSpectra[i])
						ANj = abs(st.NoiseSpectra[j])
						ANij= abs(st.NoiseSpectra[(i,j)])
						
						# Get the signal spectra for this window 
						# and subtract off known correlated noise
						y      = Y - ANij
						#~ y = it.symm_smooth( G[k], 2 ) - ANij
						y[y<0] = 0
						
						
						# this is the windowing function
						Wij = W*y/( ANi*ANj + y*(ANi+ANj) + y*y )
						# Normalize to 1, there's not a good reason 
						# to do this really.  
						if not max(Wij) == 0:	#don't divide by 0
							Wij /= max(Wij)
					elif windowing == 2:
						#PHAT
						# simply divide by the amplitude
						w = abs(G[k])
						w[w<1e-10] = 1e-10 #don't divide by 0
						Wij = W/w
					elif windowing == 3:
						#SCOT
						# uses a smooth auto correlation
						w  = it.symm_smooth( abs(ffts[i])**2, nSmooth )
						w *= it.symm_smooth( abs(ffts[j])**2, nSmooth )
						w[w<1e-10] = 1e-10 #don't divide by 0
						Wij = W/sqrt(w)
					else:
						#I don't know what this is
						print 'Unknown Windowing Option'
						sys.exit(windowing)
					
					###
					# interpolate
					g = it.fpad(G[k]*Wij, P)
					
					###
					# the cross correlation
					X =  real(fft.fftshift( fft.ifft(g) ))
					
					###
					# Normalization
					# this won't affect the peak, but does affect correlation amplitude
					# note, expensive way to calculate this
					aci = real( fft.ifft( ffts[i,:]*Wij*ffts[i,:].conj() ) )[0]
					acj = real( fft.ifft( ffts[j,:]*Wij*ffts[j,:].conj() ) )[0]
					if aci > 0 and acj > 0:
						X /= sqrt(aci*acj)
					#there's a factor of 1.5 difference between the correlation values 
					#in the 2 baseline algorithm, and here.  
					#I have not been able to track it down.  For now, so they noise 
					#reduction stays about the same for both, just correct for it.
					X  = real(X)
					###
					# Store the arrays
					xcs[k,:] = X
					bls[k]   = bl
					dls[k]   = (Settings.antDels[i] - Settings.antDels[j])/1000
					A[k,0]   = cos(th)
					A[k,1]   = sin(th)
					
					k += 1
			
			######
			# Pseudo Image
			# the pseudo image routine requires arrays as input
			###
			# compute the pseudo image
			pix   = 64
			bbox = array( [[-1.1,1.1],[-1.1,1.1]] )
			h     = cit.pseudo_image(xcs, bls, dls, A, N=pix, fs=180*P)
			
			###
			# Find the maximum pixel
			hmax  = h.argmax()/h.shape[0], h.argmax()%h.shape[1]
			
			###
			# Generate the initial solution guess
			#~ x0 = [  [(xbins[hmax[0]]+xbins[hmax[0]+1])/2],
					#~ [(xbins[hmax[1]]+xbins[hmax[1]+1])/2] ]
			x0 = [ 	[(bbox[0,1]-bbox[0,0])*(hmax[0]+.5)/pix+bbox[0,0]],
					[(bbox[1,1]-bbox[1,0])*(hmax[1]+.5)/pix+bbox[1,0]]	]
			x0 = mat(x0)	 

			######
			# Fine Timing
			# loop through the correlations and get the time of the peak 
			# accurately
			# this builds the b array
			# we'll get the xcorr amplitudes at the same time
			mMiddle = (xcs.shape[1]+1)/2
			b       = zeros( [A.shape[0], 1] )
			expeak  = zeros( A.shape[0] )
			taus = {}
			
			iR = A.shape[0]
			k = 0
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue

					#check for a record with no signal!
					if xcs[k].max() == 0:
						expeak[k] = 0
						taus[(i,j)] = 0,k
						iR -= 1
						k += 1
						continue

					###
					# Calculate time delay
					tau = (A[k,0]*x0[0,0]+A[k,1]*x0[1,0])*bls[k]/300. - dls[k]
					#convert to samples
					m0   = mMiddle - round(tau*180*P)

					###
					# parabolic fit to nearest peak
					# note, this operation is slow
					#m   = cit.find_peak(xcs[i], m)
					m   = find_peak(xcs[k], m0)

					m0  = argmax(xcs[k])
					if  (abs(m-m0) > 20*P) :#or \
						#(xcs[i,round(m)]/rms(xcs[i,:]) < 1.0) or \
						#(xcs[i,round(m)]/max(xcs[i,:]) < 0.3):
						iR -= 1

					#store the correlation value now
					###
					# there can be trouble with this if the peak finder 
					# breaks, which happens rarely
					try:
						expeak[k] = xcs[k, round(m) ]
					except:
						expeak[k] = xcs[k].max()
						iR -= 1
						#~ expeak[k] = 0
						#~ taus[(i,j)] = 0,k
						#~ iR -= 1
						#~ k += 1
						#~ continue
						

					###
					# calculate the cosine of the incident angle 
					sij  = (m/P-S)/180.
					sij -= dls[k]
					sij  = -sij
					cosi = sij*300./ bls[k]
					b[k,0] = cosi
				
					#this is for closure delay
					taus[(i,j)] = sij*1000, k
					k += 1
			
			if iR < 2:
				#~ print '454'
				#require at least 2 of the local peaks to occur reasonably 
				#close the the global peak
				self.outQ.put( [sSample, iMax, event] )
				#print 'removing - fine timing'
				continue

			######
			# Weighted Least Squares Solution
			A    = mat(A)
			b    = mat(b)
			resid = b-A*x0	

			###
			# the weighting function
			W    = eye(A.shape[0])	#the weighting function
			
			#what's the resolution of the device
			#ares = 300./80/max(bls)/2	#analog resolution
			dres = 2.2/pix/1.			#digital resolution
			#print ares, dres

			###
			# Weighted Least Squares
			continue_flag = False
			for k in range(3):
				###
				# the weighting function
				W     = eye(A.shape[0])	#the weighting function
				iR    = A.shape[0]
				eMask = ones( A.shape[0], dtype='int' )
				for i in range(A.shape[0]):
					#print abs(expeak[i]), abs(resid[i])*bls[i]/300.
					res = abs(resid[i])[0,0]	#
					if res < dres/(2*k+1)**2:
						res = dres/(2*k+1)**2
					#W[i,i] = abs(expeak[i])/abs(res*bls[i]/300.)
					W[i,i] = abs(expeak[i])**.5/abs(res)**2
					#remove solutions which are obviously bad
					if abs(b[i,0]) > 1.1:
						#print 'removing solution %i, out of bounds'%i
						#this is outside the circle
						eMask[i] = 0
						W[i,i]  *= 1e-10
						iR -= 1
					elif expeak[i] < 0:
						#print 'expeak problem'
						eMask[i] = 0
						W[i,i]  *= 1e-10
						iR -= 1
					#the typical uncertainty is around 0.01
					elif abs(resid[i,0]) > max([3*dres/(2*k+1)**2, 0.02]):
						#print '   removing solution %i, far away'%i
						#this solution is very very far away
						eMask[i] = 0
						W[i,i]  *= 1e-10
						iR -= 1
						
						
				###
				# the solution!
				# there's a chance the matrix is singular
				
				if linalg.cond(A.T*W*A) < 1./sys.float_info.epsilon:
					x1 = (A.T*W*A).I *A.T*W*(b)
					resid = b-A*x1
				else:
					#the inverse is at best tainted, not sure 
					#how this happened since the weights are limited, 
					#but it happened
					continue_flag = True
					break

			if continue_flag or iR < 3:
				#~ print '528'
				self.outQ.put( [sSample, iMax, event] )
				#~ print sSample, 'removing - baselines'
				continue
			resid = array(resid).flatten()


			######
			# Finalize the solution
			cosa,cosb = array(x1).flatten()
			
			###
			# is the solution valid? (almost imposible not to be at this point)
			if cosa**2 + cosb**2 > 1.1:
				#~ print '542'
				#do the noise stuff
				if Settings.windowing == None or Settings.noise:
					self.noiseQ.put( ffts )

				self.outQ.put( [sSample, iMax, event] )
				continue

			###
			# reduce the metrics
			# the :3 is here for compatibility reasons
			erms = mean(erms[:3] )
			epkpk = mean(epkpk[:3] )

			###
			# we should just store the resid, but doing it this way 
			# for backwards compatibility reasons
			eclos = []
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					for k in range(j+1,M):
						if j == k:
							continue
						#check that all combinations are in eMask 
						#(included in the solution)
						if 	eMask[taus[i,j][1]]==1 and \
							eMask[taus[j,k][1]]==1 and \
							eMask[taus[i,k][1]]==1:
							eclos.append( taus[i,j][0]+taus[j,k][0]-taus[i,k][0] )
			
			###
			# This filters out badly damaged solutions
			if len(eclos) == 0 or max(eclos) > 4:
				#~ print '577'
				#~ print 'closure error'
				#do the noise stuff
				if Settings.windowing == None or Settings.noise:
					self.noiseQ.put( ffts )

				self.outQ.put( [sSample, iMax, event] )
				continue
				
			eclos = min(eclos)
			#eclos  = dot( abs(resid),  bls*eMask )/.3/sum(eMask)	#the mean residual in ns
			#eclos *= 2.5
			
			###
			# The correlation peak has been problematic
			# The 0,1,2 antennas correlate better than the 3 antenna
			# so, omit correlations with this antenna
			mask = zeros( M*(M-1)/2, dtype='bool')
			k = 0
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					if i in [0,1,2] and j in [0,1,2]:
						mask[k] = True
					k += 1
			expeak= mean(expeak[mask])

			event = [erms,epkpk,expeak,eclos,cosa,cosb,iR,(iR,0,0)]
			self.outQ.put( [sSample, iMax, event] )
						
			
class WindowProcessor_2b( mp.Process):
	"""
	2 baseline processor
	"""
	def __init__( self, inQ, outQ, noiseQ, M, Settings ):
		
		#initialize mp stuff
		mp.Process.__init__(self)
		
		#initialize intf stuff
		self.inQ  = inQ		#where we get raw data
		self.outQ = outQ	#where the processed data goes
		self.noiseQ=noiseQ
		self.M    = 3		
		self.Settings = Settings
		self.running  = True
	
	def run(self):
		######
		# initialize
		Settings = self.Settings
		M       = self.M
		S  = Settings.numSamples
		I  = Settings.numIterate
		P  = Settings.numInterp
		bw = Settings.bandwidth
		windowing = Settings.windowing
		nSmooth = 5

		#calculate the time weighting, it's a modified hamming window
		#modified by being raised higher
		Wt_raise = 0.16	#for hamming, Wt_raise = 0.08
		Wt = (1-cos( 2*pi*arange(S)/S ))*(1-Wt_raise)*.5+Wt_raise

		while self.running:
			######
			# Get data from the queue
			sSample, data = self.inQ.get()

			event = None

			######
			# Generate the windowing function
			#generate the simplistic windowing function to bandwidth limit 
			#the cross correlations
			W = zeros(2*S)+1
			freq = abs( 180*fft.fftfreq( 2*S )  )
			W[ abs(freq) < bw[0] ] = 0
			W[ abs(freq) > bw[1] ] = 0


			######
			# Compute the FFT's
			ffts  = empty( (M, 2*S), dtype='complex' )
			erms  = zeros( M )
			epkpk = zeros( M )
			for i in range(M):
				d = (data[i,:] - mean(data[i,:]))*Wt
				ffts[i,:] = fft.fft( d, 2*S )
				#we track the peak of the 0 channel for overlapping windows
				if i == 0:
					iMax = argmax( abs(d) ) + sSample
				#track the amplitude
				d = real(fft.ifft( ffts[i,:]*W )[:S])
				epkpk[i] =  st.pkpk(d)
				erms[i] =  st.rms(d) 

			# Compute the crossspectra and Signal Spectra
			G   = zeros( [M*(M-1)/2, S*2], dtype='complex' )		#the cross spectra
			Wij = zeros( [M*(M-1)/2, S*2] )		#the weighting function
			Y   = zeros( [S*2] )	#the signal spectra
			k = 0
			for i in range(M):
				for j in range(i+1,M):
					if i == j:
						continue
					G [ k,:] = ffts[i,:] * (ffts[j,:].conj())
					
					if windowing in [0,2,3,None]:
						#that's none, PHAT, or SCOT
						Y[:] = 1.
					elif i in [0,1,2] and j in [0,1,2]:
						#estimate the signal spectra
						#omit antenna 3 (it's noisy)
						Y += abs(G[k])/3
					k += 1
			
			######
			# Correlate
			# We only use the 0-1 and 1-2 baselines
			A    = []	    #array, baseline angle
			xcs  = zeros( [3, S*2*P] )	#store the xcross correlations in an array
			dls  = zeros( 3 )	        #store the delays in an array
			bls  = zeros( 3 )	        #sotre the baselines in an array
			k = 0		#this is used to track the index of the arrays
			for i in range(3):
				for j in range(i+1,3):
					if i == j:
						continue
					th = it.antTheta([i,j])[0]
					bl = it.antBaseline([i,j])[0]
					
					###
					# Calculate the windowing function

					if windowing == None or windowing == 0:
						#no windowing
						Wij = W
					elif windowing == 1:
						#WII windowing
						# this one is complicated, and makes use of the
						# signal spectra Y
						ANi = abs(st.NoiseSpectra[i])
						ANj = abs(st.NoiseSpectra[j])
						ANij= abs(st.NoiseSpectra[(i,j)])
						
						# Get the signal spectra for this window 
						# and subtract off known correlated noise
						y      = Y - ANij
						y[y<0] = 0
						
						# this is the windowing function
						Wij = W*y/( ANi*ANj + y*ANi+ANj + y*y )
						# Normalize to 1, there's not a good reason 
						# to do this really.  
						if not max(Wij) == 0:	#don't divide by 0
							Wij /= max(Wij)
					elif windowing == 2:
						#PHAT
						# simply divide by the amplitude
						w = abs(G[k])
						w[w<1e-10] = 1e-10 #don't divide by 0
						Wij = W/w
					elif windowing == 3:
						#SCOT
						# uses a smooth auto correlation
						w  = it.symm_smooth( abs(ffts[i])**2, nSmooth )
						w *= it.symm_smooth( abs(ffts[j])**2, nSmooth )
						w[w<1e-10] = 1e-10 #don't divide by 0
						Wij = W/sqrt(w)
					else:
						#I don't know what this is
						print 'Unknown Windowing Option'
						sys.exit(windowing)
					
					###
					# interpolate
					g = it.fpad(G[k]*Wij, P)
					
					###
					# the cross correlation
					X =  real(fft.fftshift( fft.ifft(g) ))
					
					###
					# Normalization
					# this won't affect the peak, but does affect correlation amplitude
					# note, expensive way to calculate this
					aci = real( fft.ifft( ffts[i,:]*Wij*ffts[i,:].conj() ) )[0]
					acj = real( fft.ifft( ffts[j,:]*Wij*ffts[j,:].conj() ) )[0]
					X /= sqrt(aci*acj)
					#there's a factor of 1.5 difference between the correlation values 
					#in the 2 baseline algorithm, and here.  
					#I have not been able to track it down.  For now, so they noise 
					#reduction stays about the same for both, just correct for it.
					X  = real(X)
					###
					# Store the arrays
					xcs[k,:] = X
					bls[k]   = bl
					dls[k]   = (Settings.antDels[i] - Settings.antDels[j])/1000
					if (i == 0 and j == 1) or (i==1 and j == 2):
						A.append( [cos(th), sin(th)] )					
					k += 1
			

			######
			# Peak Finding
			# for the 2 baseline solutions, we only care where the absolute maximum is
			# also, b matrix is made to match the A matrix (only i == 0, 2)
			expeak  = zeros( 3 )
			b       = []
			taus    = zeros( 3 )
			for i in range(3):
				###
				# the guess
				m0 = argmax( xcs[i] )
				
				###
				# parabolic fit
				m = find_peak(xcs[i], m0)
				
				expeak[i] = xcs[i,m0]

				###
				# calculate the cosine of the incident angle 
				sij  = (m/P-S)/180.
				sij -= dls[i]
				sij  = -sij
				taus[i] = sij*1000	#save this for the closure calculation
				cosi = sij*300/ bls[i]
				if i == 0 or i == 2:
					b.append( [cosi] )

			######
			# Solution, by matrix match
			A    = mat(A)
			b    = mat(b)
			x    = A.I*b

			######
			# Finalize the solution
			cosa,cosb = array(x).flatten()
			###
			# is the solution valid? (almost imposible not to be at this point)
			if cosa**2 + cosb**2 > 1.1:
				#do the noise stuff
				#~ if Settings.windowing == None or Settings.noise:
					#~ st.store_spectra( ffts, ants)

				self.outQ.put( [sSample, iMax, event] )
				continue

			###
			# reduce the metrics
			# the :3 is here for compatibility reasons
			erms  = mean(erms )
			epkpk = mean(epkpk )
			expeak= mean(expeak)

			###
			# Closure delay, this is different than for 
			# the multibaseline processing
			eclos  = taus[0] + taus[2] - taus[1]
			if eclos > 4:
				self.outQ.put( [sSample, iMax, event] )
				continue				

			event = [erms,epkpk,expeak,eclos,cosa,cosb,2,(2,0,0)]
			self.outQ.put( [sSample, iMax, event] )
						
			

class WindowMerger( mp.Process):
	
	def __init__( self, outQ, Settings ):
		
		#initialize mp stuff
		mp.Process.__init__(self)
		
		#initialize intf stuff
		self.outQ     = outQ	#where the processed data goes
		self.Settings = Settings
		self.running  = True
		self.output   = []
	
	def run(self):			
		
		S = self.Settings.numSamples
		I = self.Settings.numIterate
		
		while self.running:
			self.output.append( self.outQ.get() )
			#split the output
			event, self.output = split_output( self.output, st.Settings )
			
			#sometimes there's a problem, which causes the merger to hang
			#this happens when the first event is missing a window
			#now sure how this comes up, but it does
			if len(self.output) > 2*S*I:
				self.output = self.output[1:]
			
			if event != None:
				#then we need to turn this into a real event!
				erms  = []
				epkpk = []
				expk  = []
				ecls  = []
				cosa  = []
				cosb  = []
				ants  = 0
				color = 0
				for e in event[2]:
					if e == None:
						continue
					erms. append( e[0] )
					epkpk.append( e[1] )
					expk. append( e[2] )
					ecls. append( e[3] )
					cosa. append( e[4] )
					cosb. append( e[5] )
					if e[6] > ants:
						ants = e[6]
						color=e[7]
				if ants == 0:
					continue
				if 	len(erms) == 1 and \
					not self.Settings.scribble and \
					self.Settings.numIterate < self.Settings.numSamples:
					#this is a singleton solution, they are always 
					#excluded
					continue
				eOut = st.Event( event[0], event[1], erms, epkpk, expk, ecls, 
								 cosa, cosb, ants, color )
				#print len(self.output),
				eOut.printout()
				eOut.write(self.Settings.outFileS)
				#event = st.Event(Settings.startSample,eamp0,eamp1,expeak,eclos,cosa,cosb,iR,iMax, (iR,0,0))
				#lastEvent.append( eamp0,eamp1,expeak,eclos,cosa,cosb )
				#event = [eamp0,eamp1,expeak,eclos,cosa,cosb,iR,(iR,0,0)]

#########
# Main Functions

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
		help='Path to the xns file containing the noise specta, needed for W=1 windowing')
	parser.add_argument('-T', '--threads',default=2, type=int, 
		help='The number of threads to use for processing, set this to the number of cores on the computer' )
	parser.add_argument('--scribble', default=False, action='store_true', 
		help='This makes a scribble data set, writes solution for every window' )
	parser.add_argument('--pix', default=64, type=int, 
		help='The number of pixels to use for the image' )
	parser.add_argument('-Y', '--year', default=2013, type=int,
		help='The year for the antenna arrangement' )



	#Debugging and verbose options
	parser.add_argument('-v','--verbose',action='count',
		help='Print output to std_out as well as to file')
	parser.add_argument('--noise',action='store_true',default=False,
		help='This assumes all solutions are noise, useful for getting the noise spectra only')

	helpString = """
	Examples:
	
	Low quality, quick look processing
	intf_process.py -s256 -P1 -I2048 -A3 -v -W3 [INPUT_FILE] [OUTPUT_FILE]
	
	High quality, 4 antenna processing for server grade computers
	It is not recomended that you process an entire data file in high quality
	intf_process.py -s -3770 -p 3140 -S256 -P4 -I4 -A3 -T16 -v -N [NOISE_FILE] -W1 [INPUT_FILE] [OUTPUT_FILE]
	"""

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
	st.Settings.threads   = arguments.threads
	st.Settings.scribble  = arguments.scribble
	st.Settings.pix		  = arguments.pix


	#read the input file header and store that
	st.Settings.header = st.RawHeader(st.Settings.inFileS)
	
	#from the header, we want to have the sample rate and pretrigger length on hand
	st.Settings.sampleRate  	  = st.Settings.header.sampleRate/1000000	#in samples/us
	st.Settings.preTriggerSamples = st.Settings.header.preTriggerSamples
	
	#t0 = st.Settings.header.usecond/1000
	t0 = -st.Settings.preTriggerSamples/st.Settings.sampleRate/1000
	#deal with the start and stop samples, surprising complicated, this
	if arguments.start == 0:
		#that's the beginning of the file, easy
		st.Settings.startSample = 0
	else:
		st.Settings.startSample = \
			(arguments.start-t0)*st.Settings.sampleRate*1000
	if arguments.stop == -1:	#process the whole file
		st.Settings.stopSample = -1
	else:
		st.Settings.stopSample = \
			(arguments.stop-t0)*st.Settings.sampleRate*1000
	print 'range', st.Settings.startSample, st.Settings.stopSample, arguments.start, arguments.stop
	print 't0', t0

	###
	# Deal with the year
	if arguments.year == 2012:
		# The old values used for the 2011-2012 data
		t13 = (194 - 49.08)*pi/180
		t12 = (283 - 49.08)*pi/180
		st.Settings.antLocs 	=\
			array( [ 	[ 0.00,  0.00 ],
						[ 9.88*sin(t12), 9.88*cos(t12) ],
						[10.28*sin(t13),10.28*cos(t13) ],
						[ 0.00,  0.00 ] ] )
		st.Settings.antDels	= array( [0,0,0,0] )	#equal length cables
	elif arguments.year == 2013:
		st.Settings.antLocs    =\
			array( [   	[ -5.47,-15.39 ],
						[  0.00,  0.00 ],
						[-15.97,  5.86 ],
						[-28.54,-15.63 ] ] )
		
		#Antennas delays (in ns).
		#phased together based on 2013/07/08 data
		st.Settings.antDels    =\
			array( [ 	198.00,
						210.00,
						192.03,
						162.03 ] )
	elif arguments.year == 2014:
		st.Settings.antLocs	   =\
			array( [ 	[-5.587,-14.120],
						[14.110,  3.615],
						[-8.522, 10.505],
						[ 0.000,  0.000] ] )
		st.Settings.antDels    =\
			array( [ 	198.00,
						210.00,
						192.03,
						162.03 ] )
	else:
		print 'ERROR, invalid year'
		sys.exit(arguments.year)
		
	#Are we saving or loading the noise spectra?
	if st.Settings.windowing != None and not st.Settings.noise and st.Settings.windowing != 0:
		st.Settings.storeSpec = False
	else:
		st.Settings.storeSpec = True
	
	#load the noise spectra?
	if not st.Settings.storeSpec and st.Settings.windowing not in [0,3]:
		st.read_noise_spectra()
	else:
		st.init_noise_spectra()

def split_output( output, Settings ):
	output.sort()
	
	S = Settings.numSamples
	I = Settings.numIterate
	T = Settings.threads
	
	#we need to know how many windows, at most, should contribute to an event
	Nwindows = int( S/I )
	
	#check to see if the first event is fully populated
	if len(output) < Nwindows:
		return None, output
	if output[Nwindows-1][0] - output[0][0] - S < 0:
		#this is a catch for scribble or coarse processing modes, 
		#which have only 1 window per event
		if Settings.scribble or Nwindows == 0:
			return [output[0][0], output[0][1], [output[0][2]]], output[1:]
		#there's an event to be had
		iTarget = output[0][1]
		event = [output[0][2]]
		for i in range(1,Nwindows):
			if output[i][1] - iTarget > S/2:
				break
			event.append( output[i][2] )
		
		return [output[0][0], output[0][1], event], output[i:]
	return None, output
		

def main():
	sTime = time.time()
	#get the arguments
	#build the settings
	parse_args()
	
	#~ if st.Settings.verbose:
		#~ print st.Settings.header.report()
		#~ print st.Settings.report()
	
	st.Settings.file_initialize()

	M = st.Settings.numAnts
	S = st.Settings.numSamples
	I = st.Settings.numIterate
	T = st.Settings.threads

	inQ = mp.Queue( 2*T )
	outQ = mp.Queue()
	noiseQ = mp.Queue()

	print '*** initializing threads'
	#initialize the threads
	threads = []
	for i in range( st.Settings.threads ):
		if M == 2:
			threads.append( WindowProcessor_2b( inQ, outQ,noiseQ, M, st.Settings ) )
		else:
			threads.append( WindowProcessor( inQ, outQ,noiseQ, M, st.Settings ) )
		threads[-1].start()
	
	merger = WindowMerger( outQ, st.Settings )
	merger.start()
	
	print '*** entering Main Loop'
	
	output = []
	
	#main loop
	while st.Settings.startSample <= st.Settings.stopSample or st.Settings.stopSample == -1:
		#for noise processing, this is the progress bar
		if st.Settings.noise:
			print st.Settings.stopSample-st.Settings.startSample
		
		iBuffer = 0
		#read in the data buffer
		if I < S:
			#~ print 'reading large buffer'
			dataBuffer = st.read_raw_file_data(numSamples = 16*S*T)
		else:
			#~ print 'reading small buffer'
			dataBuffer = st.read_raw_file_data(numSamples = S)
		if dataBuffer == None:
			break
		while iBuffer + S <= dataBuffer.shape[1]:
			#~ print 'appending to queue', iBuffer, dataBuffer.shape
			if inQ.full():
				time.sleep(0.001)
				continue
			###
			# Step 2: Assign the data
			#append it to the queue, this may take a while since 
			#there is a limited size to the queue
			inQ.put( [st.Settings.startSample, dataBuffer[:,iBuffer:iBuffer+S]] )

			st.Settings.startSample += I
			iBuffer += I
			if (st.Settings.startSample > st.Settings.stopSample) and \
				st.Settings.stopSample != -1:
				break

		###
		#Check the noise buffer
		while not noiseQ.empty():
			ffts = noiseQ.get()
			it.store_spectra( ffts )				

	print '*** Exiting Main Loop'

	###
	#wait for inQ to empty:
	while not inQ.empty():
		time.sleep( 0.001 )
	
	#wait for outQ to finalize
	l = outQ.qsize()
	time.sleep(0.001)
	while outQ.qsize() != l:
		l = outQ.qsize()
		time.sleep(0.001)
		print l

	print '*** Terminating Window Processors'
	#kill threads
	for t in threads:
		t.terminate()


	#wait a little longer for the merger to do it's thing
	time.sleep(.1)
	
	print '*** Terminating Window Merger'
		
	merger.terminate()	
	print '*** %i windows remain unprocessed'%len(merger.output)

	###
	#finish up the data file
	st.Settings.file_finalize()
	
	###
	# Finallay
	# store the spectrum info
	if st.Settings.windowing == None or st.Settings.noise:
		st.write_noise_spectra()
	
	print '*** completed processing in %0.2f seconds'%(time.time()-sTime)
		

if __name__ == '__main__':
	main()
