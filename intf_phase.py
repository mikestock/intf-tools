import sys,os,time,argparse
import struct

#these libraries may not exist, but are needed
try:
	from numpy import *
except:
	print 'ERROR: numpy does not appear to be installed!'
	sys.exit(2)

from scipy import optimize

import intf_tools6 as st
import intf_process6 as ip

AntOrders = [ 	[1, 0, 2],
				[3, 2, 0],
				[0, 3, 1],
				[2, 1, 3]  ]
def errfunc( p, x, y):
	return y-hyperbola(p,x)

def hyperbola( p, x ):
	#print x
	#p a,b, x0,y0
	a =p[0]
	b =p[1]
	x0=p[2]
	y0=p[3]
	#require a to be a little non-zero
	if a == 0:
		a = 1e-8
	#returns y
	return y0 + sqrt( b**2/a**2 *(x-x0)**2 +b**2 )

def main():
	st.gSettings = st.Settings( sys.argv )
	
	###The Delay
	startSample = st.gSettings.startSample
	lastEvent = None
	target = 2
	mins_t   = []
	mins_x   = []
	mins_y   = []
	

	#the delay	
	#print st.gSettings.startSample, st.gSettings.stopSample
	while st.gSettings.startSample <= st.gSettings.stopSample:

		data = st.read_ditf_file_data()
		if data == None:
			#we're out of data
			break		
		
		st.gSettings.startSample += st.gSettings.numIterate
		
		#check for saturation on all channels
		if data.max()==2**16 or data.min()==0:
			print '**** Saturation Detected ****'
			st.gSettings.SatCounter += 1
			#move forward farther than normal
			st.gSettings.startSample += st.gSettings.numSamples/2
			continue
		#check to see if it's a big source
		elif sqrt( max(data[1]**2) - min(data[1]**2) ) < 4000:
			continue

		#add a little more so the pulse is in the center of the array
		st.gSettings.startSample += 128
		
		#ok, we have a window
		
		### Delay
		fit_good = True
		phasor = []
		param = st.gSettings.antDels[target]
		params_i = arange(-20,20,2)
		for dparam in params_i:
			st.gSettings.antDels[target] = param + dparam
			lastEvent = None
			for order in AntOrders:
				event = ip.process_window(data,lastEvent, order, st.gSettings)
				if event == None:
					break
				lastEvent = event
			if lastEvent == None: continue
			solsa = array(lastEvent.cosa)
			solsb = array(lastEvent.cosb)
			cosa = sum( solsa[:,0]*solsa[:,1] )/sum(solsa[:,1])
			cosb = sum( solsb[:,0]*solsb[:,1] )/sum(solsb[:,1])
			vara = sum( solsa[:,1]*(solsa[:,0]-cosa)**2 )/sum(solsa[:,1])
			varb = sum( solsb[:,1]*(solsb[:,0]-cosb)**2 )/sum(solsb[:,1])
			estd = sqrt( vara + varb )
			phasor.append( [dparam,estd] )
		phasor = array(phasor)
		#undo the damage
		st.gSettings.antDels[target] = param

		#found all windows
		if (len(phasor) < len(params_i)):
			#this window wasn't flawless (all solutions found), keep going
			fit_good = False
		
		#min isn't on the edge
		if fit_good:
			iMin = phasor[:,1].argmin()
			if iMin == 0 or iMin == len(phasor)-1:
				#The min is at the edge of the search range, no good
				fit_good = False
		
		#least squares fit
		if fit_good:
			#do a least squares fit, 
			#initial guess
			a = phasor[iMin,1]
			b = 1.
			x0= phasor[iMin,0]
			y0= 0.
			
			fit = optimize.leastsq(errfunc,(a,b,x0,y0),args=(phasor[:,0], phasor[:,1]) )
			if fit[1] not in [1,2,3,4]:
				#the fit wasn't successful
				fit_good = False
		
		#get the quality
		if fit_good:
			#find the residual
			#y = hyperbola( (a,b,x0,y0), phasor[:,0] )
			y = hyperbola( fit[0], phasor[:,0] )			
			
			residual = sum( (phasor[:,1]-y)**2 )
			if residual > 1e-6 or abs(fit[0][3])>1e-10:
				#the fit wasn't too good
				fit_good = False
		#store result
		if fit_good:		
			qual  = fit[0][1]/fit[0][0]	#b/a, the slope of the asymptote
			delay = fit[0][2]
			
			mins_t.append( [ delay , qual] )
			print 'T',delay,x0, qual, (st.gSettings.stopSample-st.gSettings.startSample)/st.gSettings.sampleRate
		
		
		### X
		fit_good = True
		phasor = []
		param = st.gSettings.antLocs[target][0]
		params_i = arange(-5,5,1)
		for dparam in params_i:
			st.gSettings.antLocs[target][0] = param + dparam
			lastEvent = None
			for order in AntOrders:
				event = ip.process_window(data,lastEvent, order, st.gSettings)
				if event == None:
					break
				lastEvent = event
			if lastEvent == None: continue
			solsa = array(lastEvent.cosa)
			solsb = array(lastEvent.cosb)
			cosa = sum( solsa[:,0]*solsa[:,1] )/sum(solsa[:,1])
			cosb = sum( solsb[:,0]*solsb[:,1] )/sum(solsb[:,1])
			vara = sum( solsa[:,1]*(solsa[:,0]-cosa)**2 )/sum(solsa[:,1])
			varb = sum( solsb[:,1]*(solsb[:,0]-cosb)**2 )/sum(solsb[:,1])
			estd = sqrt( vara + varb )
			phasor.append( [dparam,estd] )
		phasor = array(phasor)
		#undo the damage
		st.gSettings.antLocs[target][0] = param

		#found all windows
		if (len(phasor) < len(params_i)):
			#this window wasn't flawless (all solutions found), keep going
			fit_good = False
		
		#min isn't on the edge
		if fit_good:
			iMin = phasor[:,1].argmin()
			if iMin == 0 or iMin == len(phasor)-1:
				#The min is at the edge of the search range, no good
				fit_good = False
		
		#X and Y aren't good hyperbola
		if fit_good:
			para_fit = polyfit([-1,0,1],phasor[iMin-1:iMin+2,1],2)
			para_min = [para_fit[0]*2,para_fit[1]]# the dirivative
			para_min = roots(para_min)	#the root of the derivative should be a local max/min
			if not para_min:	#the fit broke
				fit_good = False
			if abs( para_min[0] ) > 1:
				print para_min
				fit_good = False

		#store result
		if fit_good:
			qual  = phasor[:,1].max()-phasor[:,1].min()
			delay = phasor[iMin,0]+para_min[0]*(params_i[1]-params_i[0])
			
			mins_x.append( [delay , qual] )
			print 'X',delay, qual, (st.gSettings.stopSample-st.gSettings.startSample)/st.gSettings.sampleRate
		
		### Y
		fit_good = True
		phasor = []
		param = st.gSettings.antLocs[target][1]
		params_i = arange(-5,5,1)
		for dparam in params_i:
			st.gSettings.antLocs[target][1] = param + dparam
			lastEvent = None
			for order in AntOrders:
				event = ip.process_window(data,lastEvent, order, st.gSettings)
				if event == None:
					break
				lastEvent = event
			if lastEvent == None: continue
			solsa = array(lastEvent.cosa)
			solsb = array(lastEvent.cosb)
			cosa = sum( solsa[:,0]*solsa[:,1] )/sum(solsa[:,1])
			cosb = sum( solsb[:,0]*solsb[:,1] )/sum(solsb[:,1])
			vara = sum( solsa[:,1]*(solsa[:,0]-cosa)**2 )/sum(solsa[:,1])
			varb = sum( solsb[:,1]*(solsb[:,0]-cosb)**2 )/sum(solsb[:,1])
			estd = sqrt( vara + varb )
			phasor.append( [dparam,estd] )
		phasor = array(phasor)
		#undo the damage
		st.gSettings.antLocs[target][1] = param

		#found all windows
		if (len(phasor) < len(params_i)):
			#this window wasn't flawless (all solutions found), keep going
			fit_good = False
		
		#min isn't on the edge
		if fit_good:
			iMin = phasor[:,1].argmin()
			if iMin == 0 or iMin == len(phasor)-1:
				#The min is at the edge of the search range, no good
				fit_good = False
		
		#X and Y aren't good hyperbola
		if fit_good:
			para_fit = polyfit([-1,0,1],phasor[iMin-1:iMin+2,1],2)
			para_min = [para_fit[0]*2,para_fit[1]]# the dirivative
			para_min = roots(para_min)	#the root of the derivative should be a local max/min
			if not para_min:	#the fit broke
				fit_good = False
			if abs( para_min[0] ) > 1:
				print para_min
				fit_good = False

		#store result
		if fit_good:
			qual  = phasor[:,1].max()-phasor[:,1].min()
			delay = phasor[iMin,0]+para_min[0]*(params_i[1]-params_i[0])
			
			mins_y.append( [delay , qual] )
			print 'Y',delay, qual, (st.gSettings.stopSample-st.gSettings.startSample)/st.gSettings.sampleRate		
		
		#jump a lot
		st.gSettings.startSample += 12800
	
	#print the result
	mins_t = array(mins_t)
	sol_t  = sum( mins_t[:,0]*mins_t[:,1] )/sum(mins_t[:,1])
	print 'ideal delay %f based on %i windows'%(sol_t, len(mins_t))
	
	mins_x = array(mins_x)
	sol_x  = sum( mins_x[:,0]*mins_x[:,1] )/sum(mins_x[:,1])
	print 'ideal x pos %f based on %i windows'%(sol_x, len(mins_x))
	
	mins_y = array(mins_y)
	sol_y  = sum( mins_y[:,0]*mins_y[:,1] )/sum(mins_y[:,1])
	print 'ideal y pos %f based on %i windows'%(sol_y, len(mins_y))
	

	
if __name__ == '__main__':
	main()
