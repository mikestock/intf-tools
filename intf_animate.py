#!/usr/bin/python

"""intf_animate.py
produces a sequence of frames to simulate high speed video type 
imagery with the interferometer solutions
"""

########
# Written by Michael Stock, 2012
# versions: 
# the intf_tools version must match intf_process for processing to work!
#	5	-	It works!
#	6.1	-	Converted to work with intf_tool for the git repository
#			Added frame time to the title string for animations with the 
#				time frame.

#plotting libraries
from matplotlib.figure import Figure,SubplotParams
from matplotlib.font_manager import FontProperties
from matplotlib import rc, axes, gridspec
#from matplotlib.patches import Polygon
#from matplotlib.backends.backend_ps import FigureCanvasPS as FigureCanvas
#from matplotlib.backends.backend_ps import FigureCanvasPS as FigureCanvas
#from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib import cm, colors

#set dpi hight
rc('savefig',dpi=100)
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

from numpy import *
import numpy as np
import os,sys,argparse
import intf_tools as it

############
# Version
parser = argparse.ArgumentParser(description="Plot processed DITF data")
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0111')

#get the arguments
parser.add_argument('input_file',
	help='Path to ascii input file')
parser.add_argument('output_file', nargs='?', default=None,
	help='Path to png output file')
parser.add_argument('-s','--start',default=0,type=float,
	help='Start position in ms from the begininng of the file')
parser.add_argument('-p','--stop',default=-1,type=float,
	help='Stop position in ms from the beginning of the file')

#Color Options
parser.add_argument('-c','--color',  default=None,type=int,
	help='Flag for plotting time in color, instead of bw')
parser.add_argument('-a','--alpha',default=None,type=int,
	help='Which column to use for the transparancy (3,4,5)')
parser.add_argument('--alphaEx',default=0.3,type=float,
	help='Alpha Exponent for the transparancy (0-1)')

#Range Options
parser.add_argument('-A','--azimuth',default=(0,360),type=float,nargs=2,
	help='The Azimuth limits of the plot')
parser.add_argument('-E','--elevation',default=(0,90),type=float,nargs=2,
	help='The Elevation limits of the plot')
parser.add_argument('--ca',default=(-1.1,1.1),type=float,nargs=2,
	help='cos(alpha) limits' )
parser.add_argument('--cb',default=(-1.1,1.1),type=float,nargs=2,
	help='cos(beta) limits' )
parser.add_argument('-t','--time',  action='store_true',
	help='Flag for plotting time domain plots')

#Plotting Arguments
parser.add_argument('--notitles',  action='store_true',default=False,
	help='Flag for removing the titles')
parser.add_argument('--minticks',  action='store_true',default=False,
	help='Flag for minimal values along the axis')
parser.add_argument('-C','--cos',   action='store_true',
	help='Flag for plotting in the direction cosine domain' )
parser.add_argument('--animate',   default=(0.5,0.02),type=float,nargs=2,
	help='The time steps for animation (coarse, fine)' )
parser.add_argument('--intermediate',   default=[.25,0],type=float,nargs=2,
	help='[x,y] Intermediate time step x for the first y ms' )	
parser.add_argument('-N', '--nbins',   default=150,type=float,
	help='How many bins to use in the backgroudn histogram' )
parser.add_argument('-F', '--fadeframes',   default=20,type=int,
	help='Sets the persistance of frames (e-foles is this value/4)' )
parser.add_argument('-Z', '--zoomtime',   default=20,type=float,
	help='How long the zoom time window is (ms)' )
parser.add_argument('--slow',   default=None,type=str,
	help='A file of times to switch to slow plotting' )
parser.add_argument('--autoslow',  default=None, type=int,
	help='Automatically guesses when to slow down the video' )	
parser.add_argument('--autoslowtail',  default=3, type=int,
	help='How many frames to look forward for another fast thing' )	
parser.add_argument('-B','--black',   action='store_true',
	help='Plot on a Black background' )
parser.add_argument('-S', '--pointSz',   default=7.,type=float,
	help='Size of the scatter points' )
parser.add_argument('--notime', default=False, action='store_true',
	help='Plots only the left panel, no elevation vs. time')
	

	
	
#Quality Arguments
parser.add_argument('--eCls',default=1.6,type=float,
	help='Closure Error Parameter')
parser.add_argument('--eStd',default=1.6,type=float,
	help='Standard Deviation Error Parameter')
parser.add_argument('--eXpk',default=0.5,type=float,
	help='Correlation Peak Error Parameter')
parser.add_argument('--eMlt',default=0.6,type=float,
	help='Multiplicity Error Parameter')
	
parser.add_argument('-q','--qualcolumn',default=6,type=int,
	help='Which column to use for the quality (3,4,5)')

arguments = parser.parse_args(sys.argv[1:])

colorDict = {	'red': ( 	(0.0, 0.3, 0.3),
							(0.1, 0.0, 0.0),
							(0.3, 0.0, 0.0),
							(0.55, 0.0, 0.0),
							(0.8, 0.8, 0.8),
							(1.0, 1.0, 1.0) ),
				'green':(	(0.0, 0.0, 0.0),
							(0.1, 0.0, 0.0),
							(0.3, 1.0, 1.0),
							(0.55, 1.0, 1.0),
							(0.8, 1.0, 1.0),
							(1.0, 0.0, 0.0) ),
				'blue':(	(0.0, 0.5, 0.5),
							(0.1, 1.0, 1.0),
							(0.3, 1.0, 1.0),
							(0.55, 0.0, 0.0),
							(0.8, 0.0, 0.0),
							(1.0, 0.0, 0.0) ) }
#gCmap   = colors.LinearSegmentedColormap('wjet',colorDict,256)
gCmap   = it.cmap_mjet

#the gain
AntennaGain = 33	#dB
SquarePkpk  = False	#due to a bug in the processing code
minPwr      = -90	#dBm
maxPwr      = -40	#dBm


#Set the Output Datafile
if arguments.output_file == None:
	arguments.output_file = os.path.splitext(arguments.input_file)[0]+'.png'
elif os.path.splitext( arguments.output_file )[-1] == '.pdf':
	#we need a different library to write PDF's
	from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas
elif os.path.splitext( arguments.output_file )[-1] in ['.ps','.eps']:
	from matplotlib.backends.backend_ps import FigureCanvasPS as FigureCanvas
elif os.path.splitext(arguments.output_file)[1] != '.png':
	print 'ERROR: output file must be a .png'
	sys.exit(1)

#read the data
data = it.read_data_file(arguments.input_file)
if arguments.slow:
	slowTimes = loadtxt( arguments.slow, delimiter=',')
#build the times from the beginning of the second
print data.time_from_second()

data.tCls = arguments.eCls
data.tStd = arguments.eStd
data.tXpk = arguments.eXpk
data.tMlt = arguments.eMlt

data.filter()

print 'qual :',len(data.time),'solutions'


if arguments.stop != -1:
	data.tRange[1] = arguments.stop
if arguments.start != 0:
	tStart = arguments.start	#we use this to determine where to start rendering
else:
	tStart = data.tRange[0]
arguments.intermediate[1] += tStart
data.azRange = arguments.azimuth
data.elRange = arguments.elevation
data.caRange = arguments.ca
data.cbRange = arguments.cb

data.limits()

print 'limits:',len(data.time),'soloutions'




#build the primary color array
print 'color column:',arguments.color
if arguments.color == 0:
	print 'Greyscale'
	data.sort( data.time )
	#none
	c = np.zeros( (len(data.mask),4) )
elif arguments.color == 1:
	#time
	print 'By time'
	data.sort( data.time )
	c = data.time - data.time.min()
	c /= c.max()
	c = gCmap( c )
elif arguments.color == 2:
	#points
	print 'by points'
	data.sort( data.time )
	c = np.arange( len(data.mask), dtype='f' )
	c /=max(c)
	c = gCmap( c )
elif arguments.color == 3:
	#amplitude
	print 'by Amplitude'
	data.sort( data.pkpk )
	vpd = .5/2**15
	if SquarePkpk:
		intf_pkpwr = 10*log10( (data.pkpk**2/2/sqrt(2)*vpd)**2/50 )- AntennaGain +30
	else:
		intf_pkpwr = 10*log10( (data.pkpk/2/sqrt(2)*vpd)**2/50 )- AntennaGain +30
	print 'pkpk',max(intf_pkpwr), min(intf_pkpwr)
	c  = intf_pkpwr-minPwr
	c /= (maxPwr-minPwr)
	c[c<0] = 0
	c[c>1] = 1
	c = gCmap( c )

#transparency (for the primary, the transparency is reduced)
a = data.pkpk.copy()
a -= min(a)
a /= float(max(a))
print mean(a)
c[:,3] = a**arguments.alphaEx
pointSz = arguments.pointSz

#determine HistMax
aeHist = histogram2d( 	data.cosa,
						data.cosb,
						weights=data.pkpk, 
						bins=[arguments.nbins*2,arguments.nbins*2], 
						range=[data.caRange,data.cbRange] )
aeHistMax = aeHist[0].max()


#################
# Make a single plot
def bgPlot_notime(data,iT,black, tStep, bgHist, titleS=None ):
	if black:
		if tStep == arguments.animate[0]:
			txtc = 'w'
		else:
			txtc = 'gray'
		face = 'k'
		cmap = cm.bone
	else:
		if tStep == arguments.animate[0]:
			txtc = 'k'
		else:
			txtc = 'gray'
		face = 'w'	
		cmap = cm.binary

	fig = Figure( figsize=(6.4,4.8) )
	fig.subplotpars = SubplotParams(left=.10,right=.95,bottom=.125,top=.80)
	gs = gridspec.GridSpec(2, 2)
	
	
	if titleS:
		titleS += '\n%4.3f'%iT
	else:
		titleS = '%4.3f'%iT
		
	fig.text(0.5,.95, titleS, color=txtc,
			 horizontalalignment='center',
			 verticalalignment='top',)
	
	nbins = arguments.nbins

	#the primary plot, Az vs El
	mask = data.time < iT
	if arguments.cos:
		#the range rings
		theta  = linspace(0,2*pi,1000)
		circX  = sin(theta)
		circY  = cos(theta)
		aeHist = histogram2d( 	data.cosa[mask],
								data.cosb[mask],
								weights=data.pkpk[mask], 
								bins=[nbins*2,nbins*2], 
								range=[data.caRange,data.cbRange] )
		#scale the histogram
		#aeHist[0][:] -= aeHist[0].min()
		aeHist[0][:] = (aeHist[0]/aeHistMax)**.3
	else:
		aeHist = histogram2d( 	data.elev[mask],
								data.azim[mask],
								weights=data.pkpk[mask], 
								bins=[nbins*2,nbins*2], 
								range=[data.elRange,data.azRange] )
		#scale the histogram
		#aeHist[0][:] -= aeHist[0].min()
		aeHist[0][:] = (aeHist[0]/aeHistMax)**.3
	ax1 = fig.add_subplot(111,axisbg=face)

	if arguments.cos:
		ax1.plot(circX,circY,c=txtc,ls='-')
		ax1.plot(circX*cos(30*pi/180),circY*cos(30*pi/180),c=txtc,ls=':')
		ax1.plot(circX*cos(60*pi/180),circY*cos(60*pi/180),c=txtc,ls=':')
	ax1.pcolormesh( aeHist[2], aeHist[1], aeHist[0], cmap=cmap, vmax=1, vmin=0)
	ax1.yaxis.set_tick_params(labelcolor=txtc)
	ax1.yaxis.set_tick_params(color=txtc)
	ax1.xaxis.set_tick_params(labelcolor=txtc)
	ax1.xaxis.set_tick_params(color=txtc)
	if arguments.cos:
		ax1.set_xlim( data.cbRange )
		ax1.set_ylim( data.caRange )
		ax1.set_aspect('equal')
		ax1.set_xlabel('cos(alpha)', color=txtc)
		ax1.set_ylabel('cos(beta)', color=txtc)
	else:
		ax1.set_xlim( data.azRange )
		ax1.set_ylim( data.elRange )
		ax1.set_xlabel('Azimuth', color=txtc)
		ax1.set_ylabel('Elevation', color=txtc)
	ax1.locator_params(axis='x',nbins=6)

	return fig, [ax1]


def bgPlot( data,iT,black, tStep, bgHist, titleS=None ):
	if black:
		if tStep == arguments.animate[0]:
			txtc = 'w'
		else:
			txtc = 'gray'
		face = 'k'
		cmap = cm.bone
	else:
		if tStep == arguments.animate[0]:
			txtc = 'k'
		else:
			txtc = 'gray'
		face = 'w'	
		cmap = cm.binary

	fig = Figure( figsize=(6.4,4.8) )
	fig.subplotpars = SubplotParams(left=.10,right=.95,bottom=.125,top=.80)
	gs = gridspec.GridSpec(2, 2)

	if titleS:
		titleS += '\n%4.3f'%iT
	else:
		titleS = '%4.3f'%iT
	
	if titleS:
		fig.text(0.5,.95, titleS, color=txtc,
				 horizontalalignment='center',
				 verticalalignment='top',)
	
	nbins = arguments.nbins
	
	#the primary plot, Az vs El
	mask = data.time < iT
	if arguments.cos:
		#the range rings
		theta  = linspace(0,2*pi,1000)
		circX  = sin(theta)
		circY  = cos(theta)
		aeHist = histogram2d( 	data.cosa[mask],
								data.cosb[mask],
								weights=data.pkpk[mask], 
								bins=[nbins*2,nbins*2], 
								range=[data.caRange,data.cbRange] )
		#scale the histogram
		aeHist[0][:] -= aeHist[0].min()
		if aeHist[0].max() > 0:
			aeHist[0][:] = (aeHist[0]/aeHist[0].max())**0.3
	else:
		aeHist = histogram2d( 	data.elev[mask],
								data.azim[mask],
								weights=data.pkpk[mask], 
								bins=[nbins*2,nbins*2], 
								range=[data.elRange,data.azRange] )
		#scale the histogram
		aeHist[0][:] -= aeHist[0].min()
		if aeHist[0].max() > 0:
			aeHist[0][:] = (aeHist[0]/aeHist[0].max())**0.3
	ax1 = fig.add_subplot(gs[:,0],axisbg=face)
	if arguments.cos:
		ax1.plot(circX,circY,c=txtc,ls='-')
		ax1.plot(circX*cos(30*pi/180),circY*cos(30*pi/180),c=txtc,ls=':')
		ax1.plot(circX*cos(60*pi/180),circY*cos(60*pi/180),c=txtc,ls=':')
	ax1.pcolormesh( aeHist[2], aeHist[1], aeHist[0], cmap=cmap)
	ax1.yaxis.set_tick_params(labelcolor=txtc)
	ax1.yaxis.set_tick_params(color=txtc)
	ax1.xaxis.set_tick_params(labelcolor=txtc)
	ax1.xaxis.set_tick_params(color=txtc)
	if arguments.cos:
		ax1.set_xlim( data.cbRange )
		ax1.set_ylim( data.caRange )
		ax1.set_aspect('equal')
		ax1.set_xlabel('cos(alpha)', color=txtc)
		ax1.set_ylabel('cos(beta)', color=txtc)
	else:
		ax1.set_xlim( data.azRange )
		ax1.set_ylim( data.elRange )
		ax1.set_xlabel('Azimuth', color=txtc)
		ax1.set_ylabel('Elevation', color=txtc)
	ax1.locator_params(axis='x',nbins=6)

	#time axis 1
	#etHist = histogram2d( data[:,2],data[:,0],weights=data[:,4]**2, bins=[nbins/2,nbins], range=[elRange,tRange] )
	#scale the histogram
	#etHist[0][:] -= etHist[0].min()
	#etHist[0][:] = (etHist[0]/etHist[0].max())**(arguments.alphaEx/2)
	ax2 = fig.add_subplot(gs[0,1],axisbg=face)
	#ax2.pcolormesh( etHist[2], etHist[1], etHist[0], cmap=cmap)
	ax2.pcolormesh( bgHist[0][2], bgHist[0][1], bgHist[0][0], edgecolor='None',cmap=cmap)
	ax2.vlines( iT, 0, 90, 'r' )
	ax2.yaxis.set_tick_params(labelcolor=txtc)
	ax2.yaxis.set_tick_params(color=txtc)
	ax2.xaxis.set_tick_params(labelcolor=txtc)
	ax2.xaxis.set_tick_params(color=txtc)
	ax2.set_xlim( data.tRange )
	ax2.set_ylim( data.elRange )
	ax2.locator_params(axis='x',nbins=5)

	#time axis 2 (zoomed)
	#mask = data[:,0] > iT - arguments.zoomtime
	#etHist = histogram2d( data[mask,2],data[mask,0],weights=data[mask,4]**2, bins=[nbins,nbins*2], range=[elRange,[iT-arguments.zoomtime, iT]] )
	#scale the histogram
	#etHist[0][:] -= etHist[0].min()
	#etHist[0][:] = (etHist[0]/etHist[0].max())**(arguments.alphaEx/2)
	ax3 = fig.add_subplot(gs[1,1],axisbg=face)
	ax3.pcolormesh( bgHist[1][2], bgHist[1][1], bgHist[1][0],edgecolor='None', cmap=cmap)
	ax3.vlines( iT, 0, 90, 'r' )
	ax3.yaxis.set_tick_params(labelcolor=txtc)
	ax3.yaxis.set_tick_params(color=txtc)
	ax3.xaxis.set_tick_params(labelcolor=txtc)
	ax3.xaxis.set_tick_params(color=txtc)
	ax3.set_xlim( [iT-arguments.zoomtime*.75, iT+arguments.zoomtime*.25] )
	ax3.set_ylim( data.elRange )
	ax3.set_xlabel('Time (ms)', color=txtc)
	
	return fig, [ax1,ax2,ax3]


################
# Check if we're in a slow plot section
def ckSlow( iT ):
	if not arguments.slow:
		return False
	for slowTime in slowTimes:
		if iT >= slowTime[0] and iT <= slowTime[1]:
			return slowTime[0]
	return False

timeS = data.TriggerTimeS
fpsL = 1./arguments.animate[0]*1e3
fpsH = 1./arguments.animate[1]*1e3

titleS = 'Flash at %s UT\nFrame Rate: '%(timeS)

print titleS

#################
# Make the plots
iT    = tStart
frames= []
autoSlow = 0

if arguments.black:
	rc('axes', facecolor='k',edgecolor='w')

nbins = arguments.nbins
etHist = histogram2d( 	data.elev,
						data.time,
						weights=data.pkpk**2, 
						bins=[nbins/2,nbins], 
						range=[data.elRange,data.tRange] )
#Scale
etHist[0][:] -= etHist[0].min()
etHist[0][:] = (etHist[0]/etHist[0].max())**(.15)

nZoom = (data.tRange[1]-data.tRange[0])/arguments.zoomtime
etzHist = histogram2d( 	data.elev,
						data.time,
						weights=data.pkpk**2, 
						bins=[nbins/2,nbins*nZoom], 
						range=[data.elRange,data.tRange] )
#Scale
etzHist[0][:] -= etzHist[0].min()
etzHist[0][:] = (etzHist[0]/etzHist[0].max())**(0.075)

while iT < data.tRange[1]:
	if len(frames) < 10:
		iFade = iT-arguments.animate[0]*10
	else:
		iFade = frames[-10]
	#what's the time step
	#print iT, arguments.intermediate
	if iT < arguments.intermediate[1]:
		tStep = arguments.intermediate[0]
	elif ckSlow(iT+arguments.animate[0]) or ckSlow(iT):
		tStep = arguments.animate[1]
	else:
		tStep = arguments.animate[0]
	
	if arguments.autoslow:
		slowMask     = (data.time > iT) & (data.time <= iT+2*arguments.animate[0])
		slowMaskTail = (data.time > iT+2*arguments.animate[0]) & (data.time <= iT+(arguments.autoslowtail+2)*arguments.animate[0])
		#we need to catch all of the strokes
		cgMask   =  slowMask & (data.elev<10)
		#print len( data.time[slowMask] ), arguments.autoslow
		if iT < arguments.intermediate[1]:
			autoSlow = 0
		elif len(data.time[slowMask]) > arguments.autoslow:
			tStep = arguments.animate[1]
			autoSlow = 1
		elif len(data.time[cgMask]) > arguments.autoslow/4:
			tStep = arguments.animate[1]
			autoSlow = 1
		elif autoSlow > 0 and autoSlow < 3*arguments.animate[0]/arguments.animate[1]:
			tStep = arguments.animate[1]
			autoSlow += 1
		elif autoSlow > 0 and len(data.time[slowMaskTail])*3/arguments.autoslowtail > arguments.autoslow:
			print 'autoslow waiting'
			pass
		else:
			#tStep = arguments.animate[0]
			autoSlow = 0
		print len( data.time[slowMask] ), autoSlow,
	
	#build the background plot
	fpsS = '%i fps'%int( 1000./tStep )
	if arguments.notime:
		fig,ax = bgPlot_notime( data, iT, arguments.black, tStep, [etHist,etzHist], titleS+fpsS )
	else:
		fig,ax = bgPlot( data, iT, arguments.black, tStep, [etHist,etzHist], titleS+fpsS )
	
	#overlay
	mask = (data.time > iFade) & (data.time <= iT+tStep)
	iC = c[mask].copy()
	#change the opacity
	a = exp( (data.time[mask] - iT)/((iT-iFade)/4) )
	a[ a>1 ] = 1 
	iC[:,3] *= a
	
	if arguments.cos:
		ax[0].scatter( 	data.cosb[mask], 
						data.cosa[mask], 
						c=iC, marker='o',
						s=pointSz, edgecolor='None' )
	else:
		ax[0].scatter( 	data.azim[mask], 
						data.elev[mask], 
						c=iC, marker='o',
						s=pointSz, edgecolor='None' )
	
	#edit the file name
	outS = ''.join([
		os.path.splitext( arguments.output_file )[0]+'_%07i'%((iT-tStart)*1000),
		os.path.splitext( arguments.output_file )[1] ])
	print 'writing to %s'%outS
	canvas = FigureCanvas(fig)
	if arguments.black:
		canvas.print_figure(outS, dpi=150,facecolor='k',edgecolor='k')
	else:
		canvas.print_figure(outS, dpi=150)
	
	frames.append(iT)
	iT += tStep
	
	#break
	

print "To animate into an AVI, use the following command"
print "mencoder mf://%s -mf fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=1000000 -oac copy -o output.avi -ffourcc DX50"%(os.path.splitext( arguments.output_file )[0]+'_???????.png')

