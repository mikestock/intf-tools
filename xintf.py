import matplotlib
matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure, SubplotParams
from matplotlib.widgets import SpanSelector, RectangleSelector
from matplotlib import rc, gridspec, cm, colors
import numpy as np
import wx, os, gzip, time, glob,sys
from wx.lib.masked import NumCtrl


import intf_tools as it

rc('savefig',dpi=100)
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

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

class MainTab(wx.Panel):
	def __init__(self, parent,root):
		wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
		self.parent = parent
		self.root   = root

		####
		# opening buttons
		self.btnSave = wx.Button(self, label='Save')
		self.Bind(wx.EVT_BUTTON, self.OnBtnSave, self.btnSave)
		self.btnOpen = wx.Button(self, label='Open')
		self.Bind(wx.EVT_BUTTON, self.OnBtnOpen, self.btnOpen)
		self.btnNext = wx.Button(self, label='Next')
		self.Bind(wx.EVT_BUTTON, self.OnBtnNext, self.btnNext)
		self.btnPrev = wx.Button(self, label='Previous')
		self.Bind(wx.EVT_BUTTON, self.OnBtnPrev, self.btnPrev)

		self.btnReadLma = wx.Button(self, label='Read LMA')
		self.Bind(wx.EVT_BUTTON, self.OnBtnReadLma, self.btnReadLma)
		

		####
		# the only limit button, reset
		self.btnReset = wx.Button(self, label='Reset Limits')
		self.Bind(wx.EVT_BUTTON, self.OnBtnReset, self.btnReset)


		####
		# coloring buttons
		self.cmbColor = wx.ComboBox(self, 
			choices=['Greyscale','By Time','By Points','By Amplitude'],
			value='By Time', style=wx.CB_READONLY )
		self.lblColor = wx.StaticText(self, label='  Color')
		self.Bind(wx.EVT_COMBOBOX,self.OnCmbColor,self.cmbColor)
		
		self.cmbSize  = wx.ComboBox(self, 
			choices=['Small','Medium','Large','Amplitude','Exagerated'],
			value='Medium', style=wx.CB_READONLY )
		self.lblSize = wx.StaticText(self, label='  Marker Size')
		self.Bind(wx.EVT_COMBOBOX,self.OnCmbSize,self.cmbSize)

		self.cmbAlpha = wx.ComboBox(self, 
			choices=['None','Some','More'],
			value='None', style=wx.CB_READONLY )
		self.lblAlpha = wx.StaticText(self, label='  Transparency')
		self.Bind(wx.EVT_COMBOBOX,self.OnCmbAlpha,self.cmbAlpha)

		###
		# Projections
		self.cmbProjection = wx.ComboBox(self, 
			choices=['Cosine', 'Az-El'],
			value='Cosine', style=wx.CB_READONLY )
		self.lblProjection = wx.StaticText(self, label='  Projection')
		self.Bind(wx.EVT_COMBOBOX,self.OnCmbProjection,self.cmbProjection)
		
		self.topSizer = wx.BoxSizer(wx.HORIZONTAL)
		self.topSizer.Add(self.btnOpen,0,wx.RIGHT)
		self.topSizer.AddStretchSpacer(1)
		self.topSizer.Add(self.btnSave,0,wx.RIGHT)
		self.topSizer.AddStretchSpacer(1)
		self.topSizer.Add(self.btnReset,0,wx.RIGHT)
		
		
		self.btmSizer = wx.BoxSizer(wx.HORIZONTAL)
		self.btmSizer.Add(self.btnPrev,0,wx.RIGHT|wx.BOTTOM)
		self.btmSizer.AddStretchSpacer(1)
		self.btmSizer.Add(self.btnNext,0,wx.RIGHT|wx.BOTTOM)
		self.btmSizer.AddStretchSpacer(1)
		self.btmSizer.Add(self.btnReadLma,0,wx.RIGHT|wx.BOTTOM)
		
		
		
		self.grid = wx.FlexGridSizer(rows=4,cols=2,hgap=5,vgap=5)
		self.grid.Add(self.lblProjection,1,wx.LEFT)
		self.grid.Add(self.cmbProjection,1,wx.RIGHT)
		self.grid.Add(self.lblColor,1,wx.LEFT)
		self.grid.Add(self.cmbColor,1,wx.RIGHT)
		self.grid.Add(self.lblAlpha,1,wx.LEFT)
		self.grid.Add(self.cmbAlpha,1,wx.RIGHT)
		self.grid.Add(self.lblSize,1,wx.LEFT)
		self.grid.Add(self.cmbSize,1,wx.RIGHT)

		self.sizer = wx.BoxSizer(wx.VERTICAL)
		self.sizer.AddSizer(self.topSizer,0,wx.RIGHT)
		self.sizer.AddStretchSpacer(1)
		self.sizer.AddSizer(self.grid,0,wx.RIGHT)
		self.sizer.AddStretchSpacer(1)
		self.sizer.AddSizer(self.btmSizer,0,wx.RIGHT)
		
		self.SetSizer(self.sizer)
		self.Fit()

	###
	#the color comboboxes
	def OnCmbColor(self,e):
		i = e.GetSelection()
		print 'changing color option to %i'%i
		self.root.plotPanel.colorOp = i
		self.root.plotPanel.UpdatePlot()
	def OnCmbSize(self,e):
		i = e.GetSelection()
		self.root.plotPanel.sizeOp = i
		self.root.plotPanel.UpdatePlot()
	def OnCmbAlpha(self,e):
		i = e.GetSelection()
		self.root.plotPanel.alphaOp = i
		self.root.plotPanel.UpdatePlot()
	
	###
	#the projections combobox
	def OnCmbProjection(self,e):
		i = e.GetSelection()
		if i == 0:
			print 'Setting Cosine Proj.'
			self.root.plotPanel.cosine = True
		else:
			print 'Setting Az-El Proj.'
			self.root.plotPanel.cosine = False
		self.root.plotPanel.UpdatePlot()
		
	###
	#file operations
	def OnBtnReadLma(self,e):
		"""File Selector Dialog to open a data file"""
		dlg = wx.FileDialog(self, "Choose a File", "","", "*.*", wx.FD_OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			inFileS = dlg.GetPath()
		dlg.Destroy()
		self.root.OpenLma(inFileS)
		
	
	def OnBtnSave(self,e):
		defaultFileS = os.path.splitext( self.parent.parent.inFileS )[0] + '.png'
		saveFileS = None
		dlg = wx.FileDialog(self, "Choose a File", "",defaultFileS, "*.*", wx.FD_SAVE)
		if dlg.ShowModal() == wx.ID_OK:
			saveFileS = dlg.GetPath()
		dlg.Destroy()
		
		#write the file
		if saveFileS != None:
			print 'writing file',saveFileS
			self.root.plotPanel.figure_canvas.print_figure(saveFileS)
		
	def OnBtnOpen(self,e):
		"""File Selector Dialog to open a data file"""
		dlg = wx.FileDialog(self, "Choose a File", "","", "*.*", wx.FD_OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			inFileS = dlg.GetPath()
		dlg.Destroy()
		self.root.OpenFile(inFileS)

	def OnBtnNext(self,e):
		if self.root.inFileS == None:
			return
		inFileS = self.root.inFileS
		dirS = os.path.split(inFileS)[0]
		extS = os.path.splitext(inFileS)[1]
		files_i = glob.glob( '%s/LB*%s'%(dirS,extS))
		files_i.sort()
		index = files_i.index(inFileS)+1
		if index >= len(files_i):
			print 'There are no more files!!'
			return
		inFileS = files_i[index]
		self.root.OpenFile(inFileS)

	def OnBtnPrev(self,e):
		if self.root.inFileS == None:
			return
		inFileS = self.root.inFileS
		dirS = os.path.split(inFileS)[0]
		extS = os.path.splitext(inFileS)[1]
		files_i = glob.glob( '%s/LB*%s'%(dirS,extS))
		files_i.sort()
		index = files_i.index(inFileS)-1
		if index < 0:
			print 'There are no more files!!'
			return
		inFileS = files_i[index]
		self.root.OpenFile(inFileS)

	def OnBtnReset(self,e):
		self.root.plotPanel.mkPlot()

class FilterTab(wx.Panel):
	def __init__(self, parent,root):
		wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
		self.parent = parent
		self.root   = root
		
		self.font       = wx.Font(10, wx.MODERN, wx.NORMAL, wx.NORMAL, False)
		self.limitsLen  = 11
		self.filtersLen = 9
		
		###
		# Limits
		boxSize = (100,30)
		lblLimits   = wx.StaticText(self, label='  Limits   ---------------------   ')
		btnLimits   = wx.Button(self, label = 'Set Limits')
		self.Bind(wx.EVT_BUTTON, self.OnBtnLimits, btnLimits)
		
		lbltRange   = wx.StaticText(self, label='  Time')
		self.boxtRange0  = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxtRange1  = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxtRange0.SetFont(self.font)
		self.boxtRange1.SetFont(self.font)
		self.boxtRange0.fmt = '%4.3f'
		self.boxtRange1.fmt = '%4.3f'
		
		lblazRange  = wx.StaticText(self, label='  Azimuth')
		self.boxazRange0 = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxazRange1 = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxazRange0.SetFont(self.font)
		self.boxazRange1.SetFont(self.font)
		self.boxazRange0.fmt = '%3.3f'
		self.boxazRange1.fmt = '%3.3f'
		
		lblelRange  = wx.StaticText(self, label='  Elevation')
		self.boxelRange0 = wx.TextCtrl(self,wx.TE_CENTRE,size=boxSize)
		self.boxelRange1 = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxelRange0.SetFont(self.font)
		self.boxelRange1.SetFont(self.font)
		self.boxelRange0.fmt = '%2.3f'
		self.boxelRange1.fmt = '%2.3f'
		
		lblcaRange  = wx.StaticText(self, label='  cos(a)')
		self.boxcaRange0 = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxcaRange1 = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxcaRange0.SetFont(self.font)
		self.boxcaRange1.SetFont(self.font)
		self.boxcaRange0.fmt = '%1.3f'
		self.boxcaRange1.fmt = '%1.3f'
		
		lblcbRange  = wx.StaticText(self, label='  cos(b)')
		self.boxcbRange0 = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxcbRange1 = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxcbRange0.SetFont(self.font)
		self.boxcbRange1.SetFont(self.font)
		self.boxcbRange0.fmt = '%1.3f'
		self.boxcbRange1.fmt = '%1.3f'
		
		
		
		###
		# Filters
		boxSize = (80,30)
		lblFilters  = wx.StaticText(self, label='  Filters   ---------------------   ')
		btnFilters  = wx.Button(self, label = 'Set Filters')
		self.Bind(wx.EVT_BUTTON, self.OnBtnFilters, btnFilters)

		lbleCls = wx.StaticText(self, label='  eCls')
		self.boxeCls = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxeCls.SetFont(self.font)
		self.boxeCls.fmt = '%1.3f'

		lbleStd = wx.StaticText(self, label='  eStd')
		self.boxeStd = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxeStd.SetFont(self.font)
		self.boxeStd.fmt = '%1.3f'

		lbleXpk = wx.StaticText(self, label='  eXpk')
		self.boxeXpk = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxeXpk.SetFont(self.font)
		self.boxeXpk.fmt = '%1.3f'

		lbleMlt = wx.StaticText(self, label='  eMlt')
		self.boxeMlt = wx.TextCtrl(self,wx.TE_RIGHT,size=boxSize)
		self.boxeMlt.SetFont(self.font)
		self.boxeMlt.fmt = '%1.3f'
		
		
		###
		# Placement
		
		self.size1 = wx.BoxSizer(wx.HORIZONTAL)
		self.size1.Add(lblLimits,0,wx.LEFT)
		self.size1.AddStretchSpacer(1)
		self.size1.Add(btnLimits,0,wx.RIGHT)
		
		self.grid1 = wx.FlexGridSizer(rows=5,cols=3,hgap=5,vgap=5)
		self.grid1.Add(lbltRange, 1,wx.LEFT)
		self.grid1.Add(self.boxtRange0,1,wx.LEFT)
		self.grid1.Add(self.boxtRange1,1,wx.LEFT)

		self.grid1.Add(lblazRange, 1,wx.LEFT)
		self.grid1.Add(self.boxazRange0,1,wx.LEFT)
		self.grid1.Add(self.boxazRange1,1,wx.LEFT)

		self.grid1.Add(lblelRange, 1,wx.LEFT)
		self.grid1.Add(self.boxelRange0,1,wx.LEFT)
		self.grid1.Add(self.boxelRange1,1,wx.LEFT)

		self.grid1.Add(lblcaRange, 1,wx.LEFT)
		self.grid1.Add(self.boxcaRange0,1,wx.LEFT)
		self.grid1.Add(self.boxcaRange1,1,wx.LEFT)

		self.grid1.Add(lblcbRange, 1,wx.LEFT)
		self.grid1.Add(self.boxcbRange0,1,wx.LEFT)
		self.grid1.Add(self.boxcbRange1,1,wx.LEFT)

		self.size2 = wx.BoxSizer(wx.HORIZONTAL)
		self.size2.Add(lblFilters,1,wx.LEFT)
		self.size2.Add(btnFilters,0,wx.RIGHT)

		
		self.grid2 = wx.FlexGridSizer(rows=2,cols=4,hgap=5,vgap=5)
		
		self.grid2.Add(lbleCls,1,wx.LEFT)
		self.grid2.Add(self.boxeCls,1,wx.LEFT)

		self.grid2.Add(lbleStd,1,wx.LEFT)
		self.grid2.Add(self.boxeStd,1,wx.LEFT)

		self.grid2.Add(lbleXpk,1,wx.LEFT)
		self.grid2.Add(self.boxeXpk,1,wx.LEFT)

		self.grid2.Add(lbleMlt,1,wx.LEFT)
		self.grid2.Add(self.boxeMlt,1,wx.LEFT)

		self.sizer = wx.BoxSizer(wx.VERTICAL)
		self.sizer.AddSizer(self.size1,1)
		self.sizer.AddSizer(self.grid1,0)
		self.sizer.AddSizer(self.size2,1)
		self.sizer.AddSizer(self.grid2,0)
		
		self.SetSizer(self.sizer)
		self.Fit()

	###
	# Text Parser
	def text_parser(self,S):
		#first see if it's a number
		try:
			F = float(S.strip())
		except:
			return None
		return F

	###
	# Limits
	def OnBtnLimits(self,e):
		#this is ugly, I'm sure there's a better way to do it, but 
		#this should work
		data = self.root.plotPanel.data
		
		#Time
		tRange  = []
		F     = self.text_parser(self.boxtRange0.GetValue())
		if F == None:
			tRange.append( data.tRange[0] )
		else:
			tRange.append(F)
		F	  = self.text_parser(self.boxtRange1.GetValue())
		if F == None:
			tRange.append( data.tRange[1] )
		else:
			tRange.append(F)
		#Azimuth
		azRange  = []
		F     = self.text_parser(self.boxazRange0.GetValue())
		if F == None:
			azRange.append( data.azRange[0] )
		else:
			azRange.append(F)
		F	  = self.text_parser(self.boxazRange1.GetValue())
		if F == None:
			azRange.append( data.azRange[1] )
		else:
			azRange.append(F)
		#Elevation
		elRange  = []
		F     = self.text_parser(self.boxelRange0.GetValue())
		if F == None:
			elRange.append( data.elRange[0] )
		else:
			elRange.append(F)
		F	  = self.text_parser(self.boxelRange1.GetValue())
		if F == None:
			elRange.append( data.elRange[1] )
		else:
			elRange.append(F)
		#cosa
		caRange  = []
		F     = self.text_parser(self.boxcaRange0.GetValue())
		if F == None:
			caRange.append( data.caRange[0] )
		else:
			caRange.append(F)
		F	  = self.text_parser(self.boxcaRange1.GetValue())
		if F == None:
			caRange.append( data.caRange[1] )
		else:
			caRange.append(F)
		#cosb
		cbRange  = []
		F     = self.text_parser(self.boxcbRange0.GetValue())
		if F == None:
			cbRange.append( data.cbRange[0] )
		else:
			cbRange.append(F)
		F	  = self.text_parser(self.boxcbRange1.GetValue())
		if F == None:
			cbRange.append( data.cbRange[1] )
		else:
			cbRange.append(F)

		#set the limits
		self.root.plotPanel.data.tRange  =  tRange
		self.root.plotPanel.data.azRange = azRange
		self.root.plotPanel.data.elRange = elRange
		self.root.plotPanel.data.caRange = caRange
		self.root.plotPanel.data.cbRange = cbRange

		#make the changes and update the text values
		self.set_values()
		self.root.plotPanel.UpdatePlot()
	
	###
	# Filters
	def OnBtnFilters(self,e):
		F     = self.text_parser(self.boxeCls.GetValue())
		if F != None:
			self.root.plotPanel.data.tCls = F
		F     = self.text_parser(self.boxeStd.GetValue())
		if F != None:
			self.root.plotPanel.data.tStd = F
		F     = self.text_parser(self.boxeXpk.GetValue())
		if F != None:
			self.root.plotPanel.data.tXpk = F
		F     = self.text_parser(self.boxeMlt.GetValue())
		if F != None:
			self.root.plotPanel.data.tMlt = F

		#make the changes and update the text values
		self.set_values()
		self.root.plotPanel.data.filter()
		self.root.plotPanel.UpdatePlot()


	
	###
	# Reset Values
	def set_values(self):
		if self.root.plotPanel.data == None:
			return
		data = self.root.plotPanel.data
		self.boxtRange0.SetValue( (self.boxtRange0.fmt%data.tRange[0]).rjust(self.limitsLen)  )
		self.boxtRange1.SetValue( (self.boxtRange0.fmt%data.tRange[1]).rjust(self.limitsLen)  )
		self.boxazRange0.SetValue((self.boxazRange0.fmt%data.azRange[0]).rjust(self.limitsLen))
		self.boxazRange1.SetValue((self.boxazRange0.fmt%data.azRange[1]).rjust(self.limitsLen))
		self.boxelRange0.SetValue((self.boxelRange0.fmt%data.elRange[0]).rjust(self.limitsLen))
		self.boxelRange1.SetValue((self.boxelRange0.fmt%data.elRange[1]).rjust(self.limitsLen))
		self.boxcaRange0.SetValue((self.boxcaRange0.fmt%data.caRange[0]).rjust(self.limitsLen))
		self.boxcaRange1.SetValue((self.boxcaRange0.fmt%data.caRange[1]).rjust(self.limitsLen))
		self.boxcbRange0.SetValue((self.boxcbRange0.fmt%data.cbRange[0]).rjust(self.limitsLen))
		self.boxcbRange1.SetValue((self.boxcbRange0.fmt%data.cbRange[1]).rjust(self.limitsLen))

		self.boxeCls.SetValue((self.boxeCls.fmt%data.tCls).rjust(self.filtersLen))
		self.boxeStd.SetValue((self.boxeStd.fmt%data.tStd).rjust(self.filtersLen))
		self.boxeXpk.SetValue((self.boxeXpk.fmt%data.tXpk).rjust(self.filtersLen))
		self.boxeMlt.SetValue((self.boxeMlt.fmt%data.tMlt).rjust(self.filtersLen))

class CtrlPanel(wx.Notebook):
	"""Control Panel
	I know, it's a notebook
	Works largely as a container for the Tabs"""
	def __init__(self, parent,root):
		wx.Notebook.__init__(self,parent)
		#parent is gonna be important for this one
		self.parent = parent
		self.root   = root
		
		#these are the tabs
		self.fileTab = MainTab(self,self.root)
		self.AddPage(self.fileTab, "Main")

		self.filtTab = FilterTab(self,self.root)
		self.AddPage(self.filtTab, "Limit")
		
class PlotPanel(wx.Panel):
	"""Plot Panel
	Contains the plot and some plotting based functions"""
	def __init__(self, parent, root):
		wx.Panel.__init__(self,parent)
		self.parent = parent
		self.root   = root

		#important paramenters
		self.inFileS = None
		self.data    = None
		self.lma     = None
		#self.mskData = None
		self.face    = 'w'
		self.txtc    = 'k'
		self.color   = 'k'
		self.colorMap= gCmap
		self.colorOp = 1
		self.colorHL = (0,1,0,.25)
		self.alphaOp = 0
		self.sizeOp  = 1
		self.markerSz= 6
		self.marker  = 'D'
		self.cosine  = True
		self.maxA    = np.log10(65000)
		self.minA    = np.log10(  700)
		#qualities
		self.eCls = 2.25
		self.eXpk = 0.3
		self.eMlt = .40
		self.eStd = 2.0
		self._draw_pending = False
		self._draw_counter = 0
		

		self.SetBackgroundColour(wx.NamedColour("WHITE"))

		self.figure = Figure(figsize=(8.0,4.0))
		self.figure_canvas = FigureCanvas(self, -1, self.figure)

		# Note that event is a MplEvent
		self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
		#self.figure_canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)

		#this status bar is actually part of the main frame
		self.statusBar = wx.StatusBar(self.root, -1)
		self.statusBar.SetFieldsCount(1)
		self.root.SetStatusBar(self.statusBar)

		self.mkPlot()
		
		self.sizer = wx.BoxSizer(wx.VERTICAL)
		self.sizer.Add(self.figure_canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND)
		self.sizer.Add(self.statusBar, 0, wx.BOTTOM)

		self.SetSizer(self.sizer)
		self.Bind(wx.EVT_SIZE,self.OnSize)

	###
	# Makers

	def mkPlot(self):
		#clear the figure
		self.figure.clf()
		self.ax1Coll = None
		self.ax2Coll = None
		self.ax3Coll = None
		self.ax1Lma  = None
		self.ax3Lma  = None
		
		gs = gridspec.GridSpec(2, 2)

		#axis #1, the Az-El (or cosa-cosb) plot
		self.ax1 = self.figure.add_subplot(gs[:,0],axisbg=self.face)
		self.ax1.yaxis.set_tick_params(labelcolor=self.txtc)
		self.ax1.yaxis.set_tick_params(color=self.txtc)
		self.ax1.xaxis.set_tick_params(labelcolor=self.txtc)
		self.ax1.xaxis.set_tick_params(color=self.txtc)        

		#axis #2, the time-El plot (overview)
		self.ax2 = self.figure.add_subplot(gs[0,1],axisbg=self.face)
		self.ax2.yaxis.set_tick_params(labelcolor=self.txtc)
		self.ax2.yaxis.set_tick_params(color=self.txtc)
		self.ax2.xaxis.set_tick_params(labelcolor=self.txtc)
		self.ax2.xaxis.set_tick_params(color=self.txtc)        
		#self.ax2b = self.figure.add_subplot(gs[0,1],sharex=self.ax2, sharey=self.ax2, frameon=False)

		#axis #3, the time-El plot (Zoom)
		self.ax3 = self.figure.add_subplot(gs[1,1],axisbg=self.face)
		self.ax3.yaxis.set_tick_params(labelcolor=self.txtc)
		self.ax3.yaxis.set_tick_params(color=self.txtc)
		self.ax3.xaxis.set_tick_params(labelcolor=self.txtc)
		self.ax3.xaxis.set_tick_params(color=self.txtc)

		SelectorColor = self.colorHL
		self.span_ax1 = RectangleSelector(self.ax1, self.OnSelectAx1,
			minspanx=0.01, minspany=0.01,
			rectprops=dict(facecolor=self.colorHL, alpha=0.25),useblit=True)

		self.span_ax2 = SpanSelector(self.ax2, self.OnSelectAx2, 'horizontal',\
			rectprops=dict(facecolor=self.colorHL, alpha=0.25),useblit=True,minspan=0.01)

		self.span_ax3 = SpanSelector(self.ax3, self.OnSelectAx3, 'horizontal',\
			rectprops=dict(facecolor=self.colorHL, alpha=0.25),useblit=True,minspan=0.01)

		self.mkColorMap()
		

		########
		# Make the plot
		if self.data == None:
			return
		#initialize all the ranges
		self.data.reset_limits()
		self.data.update()

		print 'Making New Plot'
		
		#make the title
		self.title = self.figure.suptitle( self.data.TriggerTimeS )
		
		#Main Plot
		if self.cosine:
			print 'Cosine Projection'
			theta = np.linspace(0,2*np.pi,1000)
			X = np.cos(theta)
			Y = np.sin(theta)
			self.ax1Coll = self.ax1.scatter( 
								self.data.cosb, self.data.cosa,
								s=self.markerSz,
								marker=self.marker,
								facecolor=self.color,
								edgecolor='None' )
			self.ax1.plot(X,Y,'k-', linewidth=2)
			self.ax1.plot(np.cos(30*np.pi/180)*X,np.cos(30*np.pi/180)*Y,'k--', linewidth=2)
			self.ax1.plot(np.cos(60*np.pi/180)*X,np.cos(60*np.pi/180)*Y,'k--', linewidth=2)

			self.ax1.set_ylabel('cosb')
			self.ax1.set_xlabel('cosa')
			self.ax1.set_xlim( self.data.cbRange )
			self.ax1.set_ylim( self.data.caRange )
			self.ax1.set_aspect('equal')
		else:
			print 'Az-El Projection'
			self.ax1Coll = self.ax1.scatter( 
								self.data.azim,
								self.data.elev,
								s=self.markerSz,
								marker=self.marker,
								facecolor=self.color,
								edgecolor='None' )
			self.ax1.set_xlim( self.data.azRange )
			self.ax1.set_ylim( self.data.elRange )
			self.ax1.set_ylabel('Elevation')
			self.ax1.set_xlabel('Azimuth')
			self.ax1.set_aspect('auto')

		self.ax3Coll = self.ax3.scatter( 
							  self.data.time,
							  self.data.elev,
							  s=self.markerSz,
							  marker=self.marker,
							  facecolor=self.color,
							  edgecolor='None' )
		self.ax3.set_xlim( self.data.tRange )
		self.ax3.set_ylim( self.data.elRange )
		self.ax3.set_xlabel('Time (ms)')


		#the overview plot
		self.ax2.pcolormesh( self.data.rawDataHist[2], 
					self.data.rawDataHist[1], self.data.rawDataHist[0]**.1, 
					edgecolor='None',cmap=cm.binary)

		self.ax2Coll = 	  self.ax2.scatter( 
						  self.data.time,
						  self.data.elev,
						  s=3,
						  marker=self.marker,
						  facecolor=self.colorHL,
						  edgecolor='None' )
		#these limits shouldn't change though
		self.ax2.set_xlim( self.data.tRange  )
		self.ax2.set_ylim( self.data.elRange )

		self.root.ctrlPanel.filtTab.set_values()
		self.redraw()
		
	def mkColorMap(self):
		"""Makes a colormap"""
		print 'Color:',
		#most color maps use static sizing
		if self.data == None:
			return
		if self.colorOp == 0:
			print 'Greyscale'
			self.data.sort( self.data.time )
			#none
			self.color = np.zeros( (len(self.data.mask),4) )
		elif self.colorOp == 1:
			#time
			print 'By time'
			self.data.sort( self.data.time )
			c = self.data.time - self.data.time.min()
			c /= c.max()
			self.color = self.colorMap( c )
		elif self.colorOp == 2:
			#points
			print 'by points'
			self.data.sort( self.data.time )
			c = np.arange( len(self.data.mask), dtype='f' )
			c /=max(c)
			self.color = self.colorMap( c )
		elif self.colorOp == 3:
			#amplitude
			print 'by Amplitude'
			self.data.sort( self.data.pkpk )
			aMin = np.log10( self.data.a05 )
			aMax = np.log10( self.data.a95 )
			c = np.log10(self.data.pkpk)
			c = (c-aMin)/(aMax-aMin)
			c[c>1] = 1
			self.color = self.colorMap( c )
					
		self.mkAlpha()
	
	def mkSize(self):
		print 'MarkerSize:',
		if self.sizeOp == 0:
			#small
			print 'small'
			self.markerSz = 3
		elif self.sizeOp == 1:
			#medium
			print 'medium'
			self.markerSz = 6
		elif self.sizeOp == 2:
			#large
			print 'large'
			self.markerSz = 12
		elif self.sizeOp == 3:
			#size by amplitude
			print 'by Amplitude'
			s = np.log10( self.data.pkpk )
			s = (s-self.minA)/(self.maxA-self.minA)
			s[s>1] = 1
			s = (1+3*s**2)
			self.markerSz = 6*s
		elif self.sizeOp == 4:
			#exagerated size by ampltiude
			print 'exagerated'
			s = np.log10( self.data.pkpk )
			aMin = np.log10(self.data.aMin)
			aMax = np.log10(self.data.aMax)
			s = (s-aMin)/(aMax-aMin)
			s[s>1] = 1
			s = (1+3*s**2)**2
			self.markerSz = 6*s
			
	
	def mkAlpha(self):
		print 'Alpha:',
		if self.alphaOp == 0:
			#no alpha
			print 'None'
			self.color[:,3] = 1
			return
		elif self.alphaOp == 1:
			#some alpha
			print '0.2'
			alphaEx = .2
		elif self.alphaOp == 2:
			#more alpha
			print '0.4'
			alphaEx = .4
		else:
			#don't know this option, don't do anything
			return
		a = self.data.pkpk.copy()
		a -= min(a)
		a /= max(a)
		self.color[:,3] = a**alphaEx
		
	
	###
	#On Catches

	def OnSelectAx1(self,click, release):
		xlims = [click.xdata, release.xdata]
		xlims.sort()
		ylims = [click.ydata, release.ydata]
		ylims.sort()
		
		if self.data == None:
			return
		if self.cosine:
			self.data.caRange = ylims
			self.data.cbRange = xlims
		else:
			self.data.elRange = ylims
			self.data.azRange = xlims

		#update the plots
		self.UpdatePlot()
			
	def OnSelectAx2(self,xmin,xmax):
		self.figure_canvas.draw()
		if self.data == None:
			return
		self.data.tRange = [xmin,xmax]
		
		#update the mask and plot
		self.UpdatePlot()

	def OnSelectAx3(self,xmin,xmax):
		#mask the data array
		if self.data == None:
			return
		self.data.tRange = [xmin,xmax]
		
		#update the mask and plot
		self.UpdatePlot()

	
	def OnSize(self,e):
		if self.GetAutoLayout():
			self.Layout()
		left   = 60
		right  = 30
		top    = 30
		bottom = 40
		wspace = 100
		dpi = self.figure.dpi
		h   = self.figure.get_figheight()*dpi
		w   = self.figure.get_figwidth()*dpi
		#figure out the margins
		self.figure.subplots_adjust(left=left/w,
									right=1-right/w,
									bottom=bottom/h,
									top=1-top/h,
									wspace=wspace/w)
		self.redraw()
	
	###
	#Updaters
		
	def UpdateStatusBar(self, event):
		if event.inaxes:
			x, y = event.xdata, event.ydata
			self.statusBar.SetStatusText(( "x= " + str(x) +
										   "  y=" +str(y) ),
											0)
	#~ def UpdateMask(self):
		#~ if self.data == None:
			#~ return
		#~ 
		#~ self.data.mask = np.where( 
			#~ (self.data.time>=self.tRange[ 0])&(self.data.time<=self.tRange[ 1])&
			#~ (self.data.azim>=self.azRange[0])&(self.data.azim<=self.azRange[1])&
			#~ (self.data.elev>=self.elRange[0])&(self.data.elev<=self.elRange[1])&
			#~ (self.data.cosa>=self.caRange[0])&(self.data.cosa<=self.caRange[1])&
			#~ (self.data.cosb>=self.cbRange[0])&(self.data.cosb<=self.cbRange[1]) )[0]


	def UpdatePlot(self):
		"""redraws the main axis"""
		if self.data == None:
			return
		self.data.limits()
		self.data.update()
		self.mkColorMap()
		self.mkSize()
		
		#Main plot (remake)
		if self.ax1Coll != None:
			self.ax1Coll.remove()
		if self.ax1Lma  != None:
			self.ax1Lma.remove()
		if self.cosine:
			print 'Cosine Projection'
			self.ax1Coll = 	self.ax1.scatter( 
							self.data.cosb,
							self.data.cosa,
							s=self.markerSz,
							marker=self.marker,
							facecolor=self.color,
							edgecolor='None' )
			
			if self.lma != None:
				self.ax1Lma = self.ax1.scatter(
							self.lma.cosb,
							self.lma.cosa,
							s=6,
							marker=self.marker,
							facecolor=(0,0,0,.7),
							edgecolor='None' )
			
			self.ax1.set_ylabel('cosa')
			self.ax1.set_xlabel('cosb')
			self.ax1.set_ylim( self.data.caRange )
			self.ax1.set_xlim( self.data.cbRange )
			self.ax1.set_aspect('equal')
		else:
			print 'Az-El Projection'
			self.ax1Coll = 	self.ax1.scatter( 
							self.data.azim,
							self.data.elev,
							s=self.markerSz,
							marker=self.marker,
							facecolor=self.color,
							edgecolor='None' )
			if self.lma != None:
				self.ax1Lma = self.ax1.scatter(
							self.lma.azim,
							self.lma.elev,
							s=6,
							marker=self.marker,
							facecolor=(0,0,0,.5),
							edgecolor='None' )
			
			self.ax1.set_xlim( self.data.azRange )
			self.ax1.set_ylim( self.data.elRange )
			self.ax1.set_ylabel('Elevation')
			self.ax1.set_xlabel('Azimuth')
			self.ax1.set_aspect('auto')

		#Zoom plot (remake)
		if self.ax3Coll != None:
			self.ax3Coll.remove()
		if self.ax3Lma != None:
			self.ax3Lma.remove()
		
		self.ax3Coll = self.ax3.scatter( 
						  self.data.time,
						  self.data.elev,
						  s=self.markerSz,
						  marker=self.marker,
						  facecolor=self.color,
						  edgecolor=(1,1,1,0) )
		if self.lma != None:
			self.ax3Lma = self.ax3.scatter( 
							self.lma.time, 
							self.lma.elev, 
							s = 6,
							marker = self.marker, 
							facecolor=(0,0,0,.5),
							edgecolor='None' )
		self.ax3.set_xlim( self.data.tRange )
		self.ax3.set_ylim( self.data.elRange )
		self.ax3.set_xlabel('Time (ms)')
		
		#overview plot
		#Remake current stuff only
		if self.ax2Coll != None:
			self.ax2Coll.remove()
		
		self.ax2Coll    = self.ax2.scatter( 
						  self.data.time,
						  self.data.elev,
						  s=3,
						  marker=self.marker,
						  facecolor=self.colorHL,
						  edgecolor='None' )
		
		print "redrawing figure"
		self.root.ctrlPanel.filtTab.set_values()
		self.redraw()


	def redraw(self):
		if self._draw_pending:
			self._draw_counter += 1
			return
		def _draw():
			self.figure_canvas.draw()
			self._draw_pending = False
			if self._draw_counter > 0:
				self._draw_counter = 0
				self.redraw()
		wx.CallLater(40, _draw).Start()
		self._draw_pending = True

class MainFrame(wx.Frame):
	"""Main Frame
	All the things are stored in this"""
	def __init__(self, ):
		wx.Frame.__init__(self,None,-1, 'INTF Plot',size=(550,350))
		
		#the only paramenters stored in the frame are about the file
		self.inFileS = None
		
		#these are the main 2 panels
		self.plotPanel    = PlotPanel(self,self)	#self is the parent and root
		self.ctrlPanel    = CtrlPanel(self,self)
		
		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(self.plotPanel,1,wx.LEFT | wx.GROW)
		sizer.Add(self.ctrlPanel,0,wx.RIGHT)
		self.SetSizer(sizer)
		self.Fit()

	def OpenFile(self,inFileS):
		print 'reading data',inFileS
		self.inFileS = inFileS
		self.plotPanel.data = it.read_data_file(inFileS)
		self.plotPanel.data.time_from_second()
		print 'making plot'
		self.plotPanel.mkPlot()
		#self.plotPanel.UpdatePlot()	
	
	def OpenLma(self,inFileS):
		self.lmaFile = inFileS
		self.plotPanel.lma  = it.LmaData(inFileS)
		self.plotPanel.UpdatePlot()


class App(wx.App):

	def OnInit(self):
		'Create the main window and insert the custom frame'
		frame = MainFrame()
		self.SetTopWindow(frame)
		frame.Show(True)
		return True

if __name__=='__main__':
	app = App(0)
	app.MainLoop()
