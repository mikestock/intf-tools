import os, sys
import argparse 
import intf_tools as it
from numpy import *

#########
# VERSION
version = 'intf_xyz 5.0.130714'


useCosine = True

#constants that affect the solutions
alt 		= 3163.5	#altitude of teh interferometer

#for the first pass
aClose		= 5.0*pi/180		#how close the intf needs to be to the LMA to be used
nClose		= 20				#how many LMA sources need to be that close
nPasses     = 10				#how many passes

#for later passes
tHist		= 0.025			#how long is the history (ms)
aHist       = 0
qWeight     = -5			#negative numbers is small quality is good
minAng      = 0.2*pi/180	#minimum angle used in quality calculation
conThresh   = 50			#moves less than 100 meters in 1 pass
numThresh   = 20


def format_data( events_i, intfD ):
	lines_i = []
	for event in events_i:
		r = event[3]
		x = event[4]
		y = event[5]
		z = event[6]
		i = event[8]
		
		lineS  = '{0: >9.4f} '.format(intfD.time[i])
		# 3D location
		lineS += '{1: > 11.4f} {1: > 11.4f} {2: > 11.4f} {3: > 11.4f} '.format(r,x,y,z)
		# 2D location
		lineS += '{0: >7.2f} {1: >6.2f} {2: > 7.4f} {3: > 7.4f} '.format(intfD.azim[i], intfD.elev[i], intfD.cosa[i], intfD.cosb[i]) 
		# amplitude
		lineS += '{0: >7.1f} {1: >7.1f} '.format(intfD.pkpk[i], intfD.rms[i])
		# Sample numbers
		lineS += '{0: >9.0f} {1: >4.0f}'.format(intfD.sSam[i], intfD.iMax[i])
		#newline
		lineS += '\n'
		
		lines_i.append(lineS)
	return lines_i


#Step 1, get arguments

parser = argparse.ArgumentParser(description="Convert INTF data to 3D",fromfile_prefix_chars='@')

#file path arguments
parser.add_argument('input_intf_file',
    help='Path to the INTF file (csv)')
parser.add_argument('input_lma_file',
    help='Path to the LMA file (converted to Az/El)')
parser.add_argument('input_adx_file',nargs='?', default=None,
    help='Path to the adx file to fill in channels (optional)')
parser.add_argument('output_file', nargs='?', default=None,
    help='Path to where you want the output to go (optional)')

#flags
parser.add_argument('--version', action='version', version=version)
parser.add_argument('-v', '--verbose', action='count',
    help='Turn on Verbose output (counts)')
parser.add_argument('--debug', action='store_true',
    help='Turn on debugging output')

#paramenters
parser.add_argument('-a', '--angle', default=5.0, type=float,
    help='The searching angle, LMA sources within this distance may contribute to the range of the INTF source')
parser.add_argument('-n', '--number', default=20, type=int, 
    help='The number of LMA sources which must be withint the searching angle for the INTF source to be ranged')
parser.add_argument('-N', '--max-passes', default=10, type=int,
    help='The maximum number of passes to run')

parser.add_argument('-s', '--start', default=None, type=float,
	help='The start time of the interferometer data to be converted')
parser.add_argument('-p', '--stop', default=None, type=float, 
	help='The stop time of the interferometer data to be converted')


arguments = parser.parse_args(sys.argv[1:])

#reset some arguments based on the parser
nPasses = arguments.max_passes
aClose  = arguments.angle*pi/180
nClose  = arguments.number

intfS = arguments.input_intf_file
lmaS  = arguments.input_lma_file
print intfS
print lmaS

############
# Open the data

###
# Intf Data
intfD = it.read_data_file( intfS )
# limit to a time range
if arguments.start != None:
	intfD.tRange[0] = arguments.start
if arguments.stop != None:
	intfD.tRange[1] = arguments.stop
intfD.limits()
print 'fitting to %i sources'%(len(intfD.time))

###
# LMA Data
lmaD  = it.LmaData( lmaS )

#this is a hack!
mask =  (lmaD.cosa < 0.04) & (lmaD.z<9000) & (lmaD.cosb<0.49)
#this mask is the stuff I want to remove
mask = invert(mask)
lmaD.limits(mask)



#############
# PASS 0
print 'PASS 0'
events_i = []

# build the array of all the close events
for i in range(len(intfD.time)):
	#~ event = intD[i]
	
	#calculate the distance of the ith point to LMA points
	ang_dist = sqrt( (lmaD.cosa-intfD.cosa[i])**2+(lmaD.cosb-intfD.cosb[i])**2 )
	ang_ind = where( ang_dist < aClose )[0]

	if not len( ang_ind ) > nClose:
		#there are not enough point, continue
		continue

	r  = mean( lmaD.rang[ ang_ind ] )
	
	z = r*sin( intfD.elev[i]*pi/180 )+alt
	d = r*cos( intfD.elev[i]*pi/180 )
	x = d*sin( intfD.azim[i]*pi/180 )
	y = d*cos( intfD.azim[i]*pi/180 )
	Converged = 10000	#initial number
	Index     = i
	events_i.append( [intfD.time[i],intfD.azim[i],intfD.elev[i],r,x,y,z,Converged,Index, intfD.cosa[i], intfD.cosb[i]] )
events_i = array(events_i)

print '->', len(events_i), mean( events_i[:,7] ), max( events_i[:,7] )
outSi = 'quasi3d_pass%i.csv'%0
#build the output
lines_i = format_data( events_i, intfD )
f = open(outSi,'w')
f.writelines(lines_i)
f.close()


#############
# PASS N
iPass = 0
while iPass < nPasses:
	print 'PASS',iPass+1

	for event in events_i:
		#is it converged:
		if event[7] < conThresh/10.:
			continue

		#build a history
		history_i = events_i[ abs( events_i[:,0] - event[0] )<tHist ]

		#calculate the angular distance
		ang_dist = sqrt( (lmaD.cosa-event[9])**2+(lmaD.cosb-event[10])**2 )
		ang_ind = where( ang_dist < aClose )[0]
		
		# We calculate a quality of the solutions based on how likely 
		# the solution matches
		ang_qua = ang_dist[ ang_ind ]
		ang_qua[ ang_qua < minAng ] = minAng	#but don't weight too much
		
		#calculate the weights for the qualities based on the history
		weight = zeros(len(ang_ind))+1
		for history in history_i:
			#how far is the history from the solutions
			d = sqrt( (history[4]-lmaD.x[ang_ind])**2 +\
					  (history[5]-lmaD.y[ang_ind])**2 +\
					  (history[6]-lmaD.z[ang_ind])**2 )

			#** HACK
			#deweight solutions late in time and low in altitude
			#~ if event[0] > -1420.3:
				#~ m = lmaD.z[ang_ind] < 8800
				#~ d[m] *= 10

			weight *= (d/max(d))**(1./len(history_i))
		ang_qua *= weight
		
		#use numThresh smallest qualities
		ang_ind = ang_ind[ ang_qua.argsort()[:numThresh] ]
		ang_qua = ang_qua[ ang_qua.argsort()[:numThresh] ]
		
		#better qaulities are smaller, so qWeight is a negative exponent
		r = sum(lmaD.rang[ang_ind]*ang_qua**qWeight)/sum(ang_qua**qWeight)

		z = r*sin( event[2]*pi/180 )+alt
		d = r*cos( event[2]*pi/180 )
		x = d*sin( event[1]*pi/180 )
		y = d*cos( event[1]*pi/180 )

		#how much did this point move?
		converge = sqrt( (event[4]-x)**2 + (event[5]-y)**2 + (event[6]-z)**2 )
		
		#update the event
		event[3] = r
		event[4] = x
		event[5] = y
		event[6] = z
		event[7] = converge

	print '->', len(events_i), mean( events_i[:,7] ), max( events_i[:,7] )
	if max(events_i[:,7]) < conThresh:
		#we're done!
		break
	outSi = 'quasi3d_pass%i.csv'%(iPass+1)
	#build the output
	lines_i = format_data( events_i, intfD )
	f = open(outSi,'w')
	f.writelines(lines_i)
	f.close()
	
	iPass += 1

outSi = os.path.splitext( intfS )[0] + 'quasi3d.dat'
lines_i = format_data( events_i, intfD )
f = open(outSi,'w')
f.writelines(lines_i)
f.close()
