

import os, glob, sys, argparse
import intf_tools as it

#converts from digitizer units to V/m
Epd 	= -30./2**15
#threhold for identifying flashes
Thresh  = 1	#V/m
#slack on either side to process
Slack 	= 50	#ms

parser = argparse.ArgumentParser(description="Plot processed DITF data")
parser.add_argument('--version', action='version', version='%(prog)s 1.0.0111')

parser.add_argument('input_files', nargs='+',
	help='Path to input FA files')
	
arguments = parser.parse_args(sys.argv[1:])

arguments.input_files.sort()



for fileS in arguments.input_files:
	head = it.RawHeader( fileS )
	#the number of samples to read in
	S = (head.preTrigger+head.postTrigger)*head.samplesPerRecord
	#this will take a while
	wave = it.read_raw_waveform_file_data( fileS, head, 0, S )
	wave[1] -= 2**15
	wave[1] *= Epd
	
	m = abs(wave[1]) > Thresh
	if len( wave[0,m] ) == 0:
		continue
	t0 = round( wave[0,m][0]-Slack )
	t1 = round( wave[0,m][-1]+Slack )
	print fileS, str(t0).rjust(9), str(t1).rjust(9)
	
	
