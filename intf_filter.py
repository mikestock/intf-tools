#!/usr/bin/python

from numpy import *
import os,sys,argparse

sys.path.append( '/home/george/Documents/Dissertation/PyCodes/intf-tools')
import intf_tools as it

version = 0.1


parser = argparse.ArgumentParser(description="Removes noise from INTF data",fromfile_prefix_chars='@')
parser.add_argument('--version', action='version', version='%(prog)s'+' %f'%version)


#the implicit arguments
parser.add_argument('input_file',
	help='Path to processed input file')
parser.add_argument('output_file', nargs='?', default=None,
	help='Path to output file')

parser.add_argument('--eCls',default=1.0,type=float,
	help='Closure Delay')
parser.add_argument('--eStd',default=1.0,type=float,
	help='Standard Deviation')
parser.add_argument('--eMlt',default=0.6,type=float,
	help='Multiplicity')
parser.add_argument('--eXpk',default=0.35,type=float,
	help='Correlation Amplitude')
parser.add_argument('--elev',default=[-30,90], nargs=2,type=float,
	help='Elevation')
parser.add_argument('--azim',default=[0,360], nargs=2,type=float,
	help='Elevation')
parser.add_argument('--cosa',default=[-1.1,1.1], nargs=2,type=float,
	help='Elevation')
parser.add_argument('--cosb',default=[-1.1,1.1], nargs=2,type=float,
	help='Elevation')


arguments = parser.parse_args(sys.argv[1:])

#input output files
inFileS = arguments.input_file
if not os.path.exists(inFileS):
	print 'input file does not exist: %s'%inFileS
	sys.exit(1)
if arguments.output_file != None:
	outFileS = arguments.output_file
else:
	try:
		outFileS = inFileS.split('.dat.gz')[0]+'_filtered.dat.gz'
	except:
		outFileS = os.path.splitext(inFileS)[0] + '_filtered.dat.gz'


##############
# Read the data
intfD = it.read_data_file( inFileS )
intfD.tCls = arguments.eCls
intfD.tStd = arguments.eStd
intfD.tMlt = arguments.eMlt
intfD.tXpk = arguments.eXpk
intfD.elRange = arguments.elev
intfD.azRange = arguments.azim
intfD.caRange = arguments.cosa
intfD.cbRange = arguments.cosb
intfD.limits()
intfD.filter()
print 'Writing to %s'%outFileS
intfD.write( outFileS )
