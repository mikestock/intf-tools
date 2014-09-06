cimport numpy as np
import numpy as np
from libc.math cimport sin, cos, M_PI

def firstorder( np.ndarray[np.float64_t, ndim=1] arr, 
				np.float64_t alpha, 
				np.float64_t z=0 ):
		"""applies a first order RC filter to arr"""
		#
		
		if z==0:
			z=arr[0]
		cdef int N = arr.shape[0]
		cdef int i
		for i in range(N):
			arr[i] = alpha*arr[i]+(1-alpha)*z 
			z=arr[i]
		return arr



def find_peak( np.ndarray[ np.float64_t, ndim=1] X, int x0 ):
	#Left Peak
	cdef int iML = x0
	cdef float mML = X[x0]
	#Right Peak
	cdef int iMR = x0
	cdef float mMR = X[x0]

	cdef int i, iMax
	cdef float b,m
	cdef int RLoop2 = 1

	cdef float thresh = max(X)/10
	###
	# written for easy implementation in Cython

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
	if iML == x0:
		if X[x0] > thresh:
			RLoop2 = 0
		mML = 0
		#print 'L loop 2'
		while X[i] < thresh:
			#print i, 'guess neg'
			i -= 1
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

	while X[i] > thresh:
		#print i, 'guess pos'
		if X[i] > mMR:
			iMR = i
			mMR = X[i]
		i += 1
	#print 'R max', iMR, x0
	#Right Loop part 2
	#We don't want to do this if we started on a hill
	if iMR == x0 and RLoop2 == 1:
		mMR = 0
		#print 'R loop 2'
		while X[i] < thresh:
			i += 1
		#loop down till X goes negaitve
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
	


#~ def find_peak( np.ndarray[ np.float64_t, ndim=1] X, int x0 ):
#~ 	###
#~ 	# Fast implementation of finding and fitting a parabola to a peak
#~ 	# X is the array to find the peak of
#~ 	# x0 is an initial guess
#~ 	
#~ 	#Left Peak
#~ 	cdef int iML = x0
#~ 	cdef float mML = X[x0]
#~ 	#Right Peak
#~ 	cdef int iMR = x0
#~ 	cdef float mMR = X[x0]
#~ 
#~ 	cdef int start_on_peak = 0
#~ 	cdef int i, iMax
#~ 	cdef float b,m
#~ 
#~ 	###
#~ 	# Phase 1:
#~ 	# Choose peak candidates
#~ 
#~ 	###
#~ 	#Left Loop part 1
#~ 	i = x0-1
#~ 	while X[i] > 0:
#~ 		if X[i] > mML:
#~ 			iML = i
#~ 			mML = X[i]
#~ 		i -= 1
#~ 	#Left Loop part 2
#~ 	if iML == x0:
#~ 		start_on_peak += 1	#only one side should do this
#~ 		while X[i] < 0:
#~ 			i -= 1
#~ 		#loop down till X goes negaitve
#~ 		while X[i] > 0:
#~ 			if X[i] > mML:
#~ 				iML = i
#~ 				mML = X[i]
#~ 			i -= 1
#~ 		
#~ 	###
#~ 	#Right Loop part 1
#~ 	i = x0+1
#~ 	while X[i] > 0:
#~ 		if X[i] > mMR:
#~ 			iMR = i
#~ 			mMR = X[i]
#~ 		i += 1
#~ 	#Right Loop part 2
#~ 	if iML == x0 and start_on_peak==0:
#~ 		while X[i] < 0:
#~ 			i += 1
#~ 		#loop down till X goes negaitve
#~ 		while X[i] > 0:
#~ 			if X[i] > mMR:
#~ 				iMR = i
#~ 				mMR = X[i]
#~ 			i += 1
#~ 
#~ 	###
#~ 	# Phase 2:
#~ 	# determine which peak, right or left, is closer to the initial guess
#~ 	if x0-iML > iMR - x0:
#~ 		iMax = iML
#~ 	elif iMR - x0 > x0-iML:
#~ 		iMax = iMR
#~ 	else:	#they're equal, damnit
#~ 		if mML > mMR:
#~ 			iMax = iML
#~ 		else:
#~ 			iMax = iMR
#~ 	
#~ 	###
#~ 	# Phase 3:
#~ 	# Parabolic Fit
#~ 	#this used to be done wit built in functions, but they're apparently 
#~ 	#slow
#~ 	
#~ 	#the slopes on either side
#~ 	mML = X[iMax] - X[iMax-1]	#value taken at -1/2 (in samples)
#~ 	mMR = X[iMax+1] - X[iMax]	#value taken at +1/2 (in samples)
#~ 	
#~ 	#solving now for y = mx + b
#~ 	b = (mML + mMR)/2	#the +/- .5 will always cancel the m's
#~ 	m = mMR - mML	#the b's here alway cancel, and the .5's add to 1
#~ 	
#~ 	#and, this is the result!
#~ 	return iMax-b/m


	
	

def fourier_image( 	np.ndarray[	np.complex128_t, ndim=2] G,
					np.ndarray[	np.float64_t, ndim=1] f,
					np.ndarray[	np.float64_t, ndim=1] dl,
					np.ndarray[	np.float64_t, ndim=2] uv,
					int N=50,
					np.ndarray[ np.float64_t, ndim=2] bbox = np.array( [[-1.1,1.1],[-1.1,1.1]] ) 
					):

	cdef float l,m,theta
	cdef int h,i,j,k

	cdef np.ndarray[np.float64_t, ndim=2] Output = np.zeros( (N,N), dtype='float' )

	#calling np costs a lot, so lets get that out of the loop as much as possibe
	#cdef np.ndarray[np.float64_t, ndim=2] Gr = np.real(G)
	#cdef np.ndarray[np.float64_t, ndim=2] Gi = np.imag(G)
	#cdef float pi = np.pi
	
	#Loop through pixels
	for i in range(N):
		l = (bbox[0,1]-bbox[1,0])*(i+.5)/N+bbox[0,0]
		for j in range(N):
			m = (bbox[1,1]-bbox[0,0])*(j+.5)/N+bbox[1,0]
#~ 			if (m**2+l**2) > 1.1:
#~ 				#over the horizon
#~ 				continue			
			#loop through baselines
			for k in range( G.shape[0] ):
				
				#loop through frequency
				#we're running out of indices!
				for h in range(f.shape[0]):
					if f[h] < 0:
						continue
					theta = -2*M_PI*f[h]*(uv[k,0]*l+uv[k,1]*m-dl[k])
					
					Output[i,j] += 2*G[k,h].real*cos(theta)
					Output[i,j] -= 2*G[k,h].imag*sin(theta) 
					
	return Output

def proj_image( np.ndarray[	np.float64_t, ndim=2] xc, 
					np.ndarray[	np.float64_t, ndim=1] bl, 
					np.ndarray[	np.float64_t, ndim=1] dl,
					np.ndarray[ np.float64_t, ndim=2] A, 
					int N=50, 
					float fs=360,
					np.ndarray[ np.float64_t, ndim=2] bbox = np.array( [[-1.1,1.1],[-1.1,1.1]] )
					):
	

	cdef float cosa, cosb, tau, dtau

	#these are used to calculate the amplitude of the pixel
	cdef float p, l
	
	#counters, there's a bunch of them
	cdef int i,j,k,m

	#int, corresponds to tau=0
	cdef int  mMiddle = (xc.shape[1]+1)/2	

	#the size of the grid
	cdef float dcosa = (bbox[0,1]-bbox[0,0])/N*.5
	cdef float dcosb = (bbox[1,1]-bbox[1,0])/N*.5

	###
	# Generate the image array
	# can't be empty since not all values will be calculated
	cdef np.ndarray[np.float64_t, ndim=2] Output = np.zeros( (N,N), dtype='float' )
	
	###
	# We're gonna do this the slow way, which is gonna take a while
	# Loop through the pixel
	for i in range(N):
		cosa = (bbox[0,1]-bbox[0,0])*(i+.5)/N+bbox[0,0]	
		for j in range(N):	
			cosb = (bbox[1,1]-bbox[1,0])*(j+.5)/N+bbox[1,0]
			
			###
			# loop over the baselines
			for k in range( xc.shape[0] ):
				tau = (A[k,0]*cosa+A[k,1]*cosb)*bl[k]/300. - dl[k]
				tau = tau*fs

				dtau = fs*dcosa*bl[k]/300.				
				
				if dtau > .5:

					m = int(tau-dtau)	#this rounds down
					
					#there's a bunch of tau values that will contibute
					l = 0
					p = 0
					while m < tau+dtau:
						#there's a chance that we're out of range on the xc
						if mMiddle-m < 0 or mMiddle-m >= xc.shape[1]:
#~ 							print 'out of bounds', k, mMiddle-m
							m += 1
							continue
						if m < tau-dtau:
#~ 							print 'l',(1-tau+dtau+m)
							l += (1-tau+dtau+m)
							p += xc[k,mMiddle-m]*(1-tau+dtau+m)
						elif m+1 > tau+dtau:
#~ 							print 'u',(tau+dtau-m)
							l += (tau+dtau-m)
							p += xc[k,mMiddle-m]*(tau+dtau-m)
						else:
							l += 1
							p += xc[k,mMiddle-m]						
						m += 1
					if l != 0:
						Output[i,j] += p/l
				else:
					m = int(tau)
					if mMiddle-m-1 < 0 or mMiddle-m >= xc.shape[1]:
						continue
					#there are just edge values, use linear interpolation
					l = tau - m

					#the 2 additions take into account some linear interpolation
					Output[i,j] += (1-l)*xc[k,mMiddle-m]
					Output[i,j] +=    l *xc[k,mMiddle-m-1]				

	return Output

#being kept for legacy reasons
def pseudo_image( 	np.ndarray[	np.float64_t, ndim=2] xc, 
					np.ndarray[	np.float64_t, ndim=1] bl, 
					np.ndarray[	np.float64_t, ndim=1] dl,
					np.ndarray[ np.float64_t, ndim=2] A, 
					int N=50, 
					float fs=360,
					np.ndarray[ np.float64_t, ndim=2] bbox = np.array( [[-1.1,1.1],[-1.1,1.1]] )
					):	
	
	cdef float cosa, cosb, tau
	cdef int i,j,k,m
	
	cdef int  mMiddle = (xc.shape[1]+1)/2	#int, corresponds to tau=0
	
	###
	# Generate the image array
	# can't be empty since not all values will be calculated
	cdef np.ndarray[np.float64_t, ndim=2] Output = np.zeros( (N,N), dtype='float' )
	
	###
	# We're gonna do this the slow way, which is gonna take a while
	# Loop through the pixel
	#with nogil:
	for i in range(N):
		for j in range(N):
			#cosa = 2*(i-N/2.+arange(2))/N
			#cosb = 2*(j-N/2.+arange(2))/N
			#cosa = 2.2*(i-N/2+.5)/N
			#cosb = 2.2*(j-N/2+.5)/N
			
			cosa = (bbox[0,1]-bbox[0,0])*(i+.5)/N+bbox[0,0]
			cosb = (bbox[1,1]-bbox[1,0])*(j+.5)/N+bbox[1,0]
			
#~ 			if (cosa**2+cosb**2) > 1.1:
#~ 				#over the horizon
#~ 				continue
			#loop through the xcorrs
			for k in range( xc.shape[0] ):
				###
				# determine how many samples off middle we are
				tau = (A[k,0]*cosa+A[k,1]*cosb)*bl[k]/300. - dl[k]
				tau = tau*fs
				
				m = int(tau)	#this rounds down
				
				#print m, tau
				
				#now I need the fractional part of tau
				tau = tau - m
				
				#there's a chance that we're out of range on the xc
				if mMiddle-m-1 < 0 or mMiddle-m >= xc.shape[1]:
					continue
				
				#the 2 additions take into account some linear interpolation
				Output[i,j] += (1-tau)*xc[k,mMiddle-m]
				Output[i,j] +=    tau *xc[k,mMiddle-m-1]
#~ 				m = int( round(tau) )
#~ 				Output[i,j] +=    xc[k,mMiddle-m]
				

	return Output
	
def mult_base_gaussian( np.ndarray[	np.float64_t, ndim=2] h, 
						np.ndarray[	np.float64_t, ndim=1]  x0, 
						float sigma, 
						np.ndarray[	np.float64_t, ndim=1] xbins, 
						np.ndarray[	np.float64_t, ndim=1] ybins ):
	
	cdef int i, j
	cdef float r
	for i in range(h.shape[0]):
		for j in range(h.shape[1]):
			r = ((xbins[i]+xbins[i+1])/2 - x0[0])**2 + ((ybins[j]+ybins[j+1])/2 -x0[1])**2
			h[i,j] *= np.exp( -r/2/sigma**2 )
	
	return h
