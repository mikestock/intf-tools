intf-tools
==========

Tools for processing, displaying and simulating lightning interferometer data
These tools were initially developed for the DITF developed at New Mexico 
Inst. of Mining and Technology as part of a PhD thesis by Michael Stock.  

The tools:

intf_tools - 
	This is the main library of the package.  It contains functions for 
	reading raw and processes data, as well as tools for simulating raw 
	waveforms and filting signal.  

intf_process - 
	Takes raw data (voltage waveforms) and makes it processed data 
	(azimuth, elevation, time).  There are many command line options 
	which.

intf-animate -
	Takes processed data and makes a sequence of frames to generate a movie

xintf	-	
	wxPython GUI.  The plotting backend is still matplotlib, which means 
	that scatter is rather slow and clunky.  But it does work.

intf_phase - 
	This is a hacked together tool which needs revision.  It is there 
	to help determine the proper locations and delays of the antennas 
	by making processing data and having the sources converge.  The 
	process is currently not automatic at all.
