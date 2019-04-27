import numpy as np
import re
from __casac__ import *
from taskinit import casalog
import delmod
import calRcore.robustcal
rantsol=calRcore.robustcal.rantsol
import time as clock

import sys
import VisAccum

class calibSolver:
	def __init__(self,vis,caltable,field='',spw='',uvrange='',scan='',observation='',
	             solint='inf',combine='',
	             refant=0,minblperant=6,minsnr=5,calmode='ap',outlieralgorithm='RMS',
                     threshold=[7,5,3],dosolintflag=False,mingain=1e-2,
		     gaintable=[],gainfield=[],interp=[],spwmap=[],debugmode=True):
		  
		self.debug=debugmode
		self.ms=ms.ms()
		self.ms2=ms.ms()
		self.tb=table.table()		
		self.vis=vis
		self.Nant=self.getNant()
		self.Nb=(self.Nant*(self.Nant-1))//2
		self.solint=self.parseSolInt(solint)
		self.field=field
		self.scan=scan
		self.uvrange=uvrange
		self.calmode=calmode
		self.caltable=caltable
		self.spw=spw
		self.observation=observation
		self.minsnr=float(minsnr)
		self.startrow=0
		self.threshold=threshold
		self.uvrange=uvrange
		self.refant=refant
		self.minblperant=minblperant
		self.outlieralgorithm=outlieralgorithm
		self.dosolintflag=dosolintflag
		self.combine=combine
		self.error=False
		self.mingain=mingain
		self.gaintable=gaintable
		self.gainfield=gainfield
		self.interp=interp
		self.spwmap=spwmap
		self.modelrestore=False
		self.dcolumn='DATA'
		self.tblock=300.0
		self.ispw=0	
		self.error=False	
		self.nCorr,self.totalIpoln=self.getPolnTypes()	
		casalog.post("CalR_v2.5")
		if(self.debug):
			self.initDebugCounters()			
		else:
			np.seterr(divide='ignore', invalid='ignore')
			np.warnings.filterwarnings('ignore')
		#casalog.origin('gaincalR')
		#self.createCaltable()
		if(self.solint !='inf' and self.solint != 'int'):
			self.solint=self.solint*60.0
		self.dataDescIndex,self.spwIndex=self.getSpwMap() 
		self.solintMap=self.getSolintMap() 
		self.initVisAccum()
		#solver statistics counters
		self.nflaginit=0
		self.nflagsolver=0
		self.ntotalsol=0
		self.nflagfinal=0
		

		
	def initDebugCounters(self):
		self.timemsio=0.0
		self.timeaccum=0.0
		self.timemean=0.0
		self.timesolve=0.0
	
	def createCaltable(self,tableType,singleChan):
		cb=calibrater.calibrater()
		cb.open(self.vis,False,False,False)
		cb.selectvis(field=self.field,spw=self.spw)
		cb.createcaltable(self.caltable, 'Complex', tableType, singleChan) 
		cb.close()
	
	def getNant(self):
		#Not reading from ANTENNA table because it contains 2 extra ant for GMRT
		#This will mess up the assumptions I make about data sorting.
		self.tb.open(self.vis+'/ANTENNA',nomodify=True)
		
		nant=self.tb.nrows()
		self.tb.close()
		if(self.debug):
			casalog.post("DEBUG: Number of antennas: %d"%nant)
		return nant
		
	def parseSolInt(self,solint):
		error=False
		if(solint!='int' and solint!='inf'):
			matchObj=re.match("(.*)s",solint)			
			if(matchObj!=None and len(matchObj.groups())==1):
				try:
				    solint=float(matchObj.group(1))
				    solint/=60.0
				except:
				    error=True		
				
			else:
				matchObj=re.match("(.*)min",solint)
				if(matchObj!=None and len(matchObj.groups())==1):
					try:
					    solint=float(matchObj.group(1))
					except:
					    error=True			
				else:
					error=True
			if(not error and self.debug):
				casalog.post("DEBUG: solution interval = %.6f min"%solint)
				
		if(error==True):
			self.error=(error or self.error)
			casalog.post("Cannot understand solint. E.g. 10s, 2min, inf, int")	
			return -1

			
		return solint
	
	 	
	
	def corrupt(self,scanrange):
		if(len(self.gaintable)==0):
			return
		casalog.post(scanrange)
    		casalog.post('Backing up MODEL_DATA column')
    		if(not self.hasmodel):
    			self.hasmodel=True
    			self.keystoget.append('MODEL_DATA')     
		self.modelrestore=True
		self.ms2.open(self.vis,nomodify=True)
		self.ms2.msselect({'SCAN':scanrange,'FIELD':self.field})
		if(self.hasmodelinit):
			self.msbackup=self.ms2.getdata(['MODEL_DATA','FLAG_ROW','FLAG'])
		else:
			self.msbackup=self.ms2.getdata(['FLAG_ROW','FLAG'])
		self.ms2.close()
							
		cb=calibrater.calibrater()
		cb.open(self.vis)
		cb.selectvis(scan=scanrange,field=self.field,spw=self.spw)
		#The following part was copy pasted from the casa gaincal task (see task_gaincal.py)
		
		ngaintab = 0;
		if (self.gaintable!=['']):
			ngaintab=len(self.gaintable)

		ngainfld = len(self.gainfield)
		nspwmap = len(self.spwmap)
		ninterp = len(self.interp)

		# handle list of list issues with spwmap
		if (nspwmap>0):
			if (type(self.spwmap[0])!=list):
				# first element not a list, only one spwmap specified
				# make it a list of list
				self.spwmap=[self.spwmap];
				nspwmap=1;

		for igt in range(ngaintab):
			if (self.gaintable[igt]!=''):

				# field selection is null unless specified
				thisgainfield=''
				if (igt<ngainfld):
					thisgainfield=self.gainfield[igt]
				
				# spwmap is null unless specifed
				thisspwmap=[-1]
				if (igt<nspwmap):
					thisspwmap=self.spwmap[igt];

				# interp is 'linear' unless specified
				thisinterp='linear'
				if (igt<ninterp):
					if (interp[igt]==''):
						interp[igt]=thisinterp;
					thisinterp=self.interp[igt];
				#calwt False is the only modification from CASA's gaincal code
				#This is because of the uncertainty with how weights are used
				#in gaincal at this moment.
				cb.setapply(t=0.0,table=self.gaintable[igt],field=thisgainfield,
					      calwt=False,spwmap=thisspwmap,interp=thisinterp)
		cb.corrupt()
		cb.correct('flag')
		cb.done()
		cb.close()
	
	def restoremodel(self,scanrange):
		
		if(not self.modelrestore):
			return
		casalog.post('Restoring MODEL_DATA column')
		self.modelrestore=False
		self.ms2.open(self.vis,nomodify=False)
		self.ms2.msselect({'SCAN':scanrange,'FIELD':self.field})
		self.ms2.putdata(self.msbackup)
		self.ms2.close()
		
	
	def getNant(self):
		#Not reading from ANTENNA table because it contains 2 extra ant for GMRT
		#This will mess up the assumptions I make about data sorting.
		self.tb.open(self.vis+'/ANTENNA',nomodify=True)
		
		nant=self.tb.nrows()
		self.tb.close()
		if(self.debug):
			casalog.post("DEBUG: Number of antennas: %d"%nant)
			
		self.ant1=np.zeros((nant*(nant-1)/2),dtype=np.int32)
		self.ant2=np.zeros((nant*(nant-1)/2),dtype=np.int32)
		ibl=0
		for iant1 in range(0,nant):
			for iant2 in range(iant1+1,nant):
				self.ant1[ibl]=iant1
				self.ant2[ibl]=iant2
				ibl+=1
		return nant
		
		
		
	def putInTable(self,gains,gains_er,snr,flags,time,scan,field,spw):
			
		nrows=flags.shape[2]
		one=np.ones(nrows)
		self.tb.open(self.caltable,nomodify=False)
		self.tb.addrows(nrows)
		self.tb.putcol('TIME',one*time,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('FIELD_ID',one*field,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('SPECTRAL_WINDOW_ID',one*spw,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('ANTENNA1',np.arange(0,nrows),startrow=self.startrow,nrow=nrows)
		self.tb.putcol('ANTENNA2',one*-1,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('SCAN_NUMBER',one*scan,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('OBSERVATION_ID',one*0,startrow=self.startrow,nrow=nrows) #Forcing this to zero. Cannot retrive obs_id. May cause issues later.
		self.tb.putcol('CPARAM',gains,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('PARAMERR',gains_er,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('FLAG',flags,startrow=self.startrow,nrow=nrows)
		self.tb.putcol('SNR',snr,startrow=self.startrow,nrow=nrows)
		self.tb.flush()
		self.tb.close()
		self.startrow+=nrows
	
	def getInc(self,time,scan,index):
		if(self.solint=='int'):
			return 1
		elif(self.combine==''):
			if(self.solint=='inf'):
				return (np.argwhere(scan[index:]==scan[index])[-1][0]+1)
			else:
				if(index+self.maxSamples>=len(scan)):
					return (len(scan)-index)
				if(scan[index]==scan[index+self.maxSamples]):
					return self.solintSamples
				else:
					return (np.argwhere(scan[index:]==scan[index])[-1][0]+1)
			
		elif(self.combine=='scan'):
			if(self.solint=='inf'):
				return scan[-1]
			else:
				return np.argwhere(time[index:]>time[index]+self.solint)[-1][0]+1
	
	def getPolnTypes(self):
		self.tb.open(self.vis+"/DATA_DESCRIPTION")
		polSetup=np.unique(self.tb.getcol("POLARIZATION_ID"))
		self.tb.close()
		self.tb.open(self.vis+"/POLARIZATION",nomodify=True)
		for iPolSetup in polSetup:		
	      		#self.ms.selectinit()  
			corr_types=self.tb.getcell('CORR_TYPE',rownr=iPolSetup)
			poln=''
			nCorr=0
			for corr_type in corr_types:
				if(corr_type==5):
					poln+='RR,'
					nCorr+=1
				elif(corr_type==8):
					poln+='LL,'
					nCorr+=1
				elif(corr_type==9):
					poln+='XX,'
					nCorr+=1
				elif(corr_type==12):
					poln+='YY,'
					nCorr+=1
		self.tb.close()
		poln=poln[:-1]
		if(self.debug):
			casalog.post("Selecting polarizations: %s"%poln)
		return nCorr,poln
		
	def getSpwMap(self):
		self.tb.open(self.vis+"/DATA_DESCRIPTION",nomodify=True)
      		#self.ms.selectinit()  
		spwMS=self.tb.getcol('SPECTRAL_WINDOW_ID')
		if(self.debug):
			casalog.post("Spectral windows in ms:")
			casalog.post(np.array2string(spwMS))
		self.tb.close()
		dataDescIndex=[]
		spwIndex=[]
		for i in range(0,len(spwMS)):
			self.ms.open(self.vis,nomodify=True)
			self.ms.selectinit(datadescid=i)
			self.ms.msselect({'field':self.field,'spw':self.spw,'scan':self.scan,'uvdist':self.uvrange,
				'observation':self.observation})			
			if(self.ms.nrow(selected=True)):
				dataDescIndex.append(i)
				spwIndex.append(spwMS[i])
			self.ms.close()
		if(self.debug):
			casalog.post("Selected DATA_DESCRIPTION_ID:")
			casalog.post(np.array2string(np.array(dataDescIndex)))
		return dataDescIndex,spwIndex
		
	def getSolintMap(self):
		solintMap=[]
		for ispw in range(len(self.spwIndex)):
			thisSolintMap=[]
			self.ms.open(self.vis,nomodify=True)
	      		self.ms.selectinit(datadescid=self.dataDescIndex[ispw])  
			self.ms.msselect({'field':self.field,'spw':self.spw,'scan':self.scan,'uvdist':self.uvrange,
					'observation':self.observation})
			dataAll=self.ms.getdata(['TIME','SCAN_NUMBER'],ifraxis=True)
			scan=dataAll['scan_number']	
			time=dataAll['time']	
			if(self.solint!='inf' and self.solint!='int'):		
				self.solintSamples=int((self.solint+0.1)/(time[1]-time[0]))
				self.maxSamples=int((self.solint*1.6+0.1)/(time[1]-time[0]))+1

			i=0
			while(i<scan.shape[0]):
				inc=self.getInc(time,scan,i)
				thisSolintMap.append(inc)
				i+=inc
				#casalog.post(i,self.solintSamples,self.maxSamples,inc)
			self.ms.close()
			solintMap.append(thisSolintMap)
		solintMap=np.array(solintMap)
		if(self.debug):
			casalog.post("DEBUG: Solint map")
			casalog.post(np.array2string(solintMap))
		return solintMap
			
	
	def initVisAccum(self):
		self.accum=[]
		for ispw in range(len(self.dataDescIndex)):
			#if(self.debug):
			#	casalog.post(self.vis,self.field,self.spw,self.uvrange,self.scan,self.observation,np.max(self.solintMap[ispw]),float(self.tblock),int(self.dataDescIndex[ispw]),int(self.spwIndex[ispw]))		
			self.accum.append(VisAccum.VisAccum(self.vis, self.field,self.spw, self.uvrange,self.scan,self.observation,self.totalIpoln,np.max(self.solintMap[ispw]), float(self.tblock),int(self.dataDescIndex[ispw]),int(self.spwIndex[ispw])))
		
		

			
	def attachMS(self):
		self.ms.open(self.vis,nomodify=True)
      		self.ms.selectinit(datadescid=self.dataDescIndex[self.ispw])  
		self.ms.msselect({'field':self.field,'scan':self.scan,'uvdist':self.uvrange,
				'observation':self.observation,'spw':'%d'%self.spwIndex[self.ispw],'correlation:':self.totalIpoln})
		#if(len(self.chanindex)!=0):	
		#	self.ms.selectchannel(start=self.chanindex[0,1],nchan=self.chanindex[0,2]-self.chanindex[0,1],inc=self.chanindex[0,3])
				
	def coreSolve(self,accumData,accumModel,accumWeight,accumFlag):
		ispw=self.ispw
		accumFlag=np.logical_not(accumFlag)
		accumWeight=1.0/accumWeight
		goodbl=(np.mean(np.mean(np.mean(accumFlag,axis=0),axis=1),axis=1)>0.1)
		casalog.origin("gaincalRIO::coreIO::coreSolve")		
		if(self.debug):
			casalog.post("Good baselines: %d"%np.sum(goodbl))
		solver=rantsol(threshold=self.threshold,minblperant=self.minblperant,refant=self.refant,
		     			      Nant=self.Nant,ant1=self.ant1[goodbl],ant2=self.ant2[goodbl],
		       			      outlieralgorithm=self.outlieralgorithm,
		       			      dorobustmean=False,calmode=self.calmode,mingain=self.mingain,
		       			      minsnr=self.minsnr,dosolintflag=self.dosolintflag, 						      						      debug=self.debug)
				
		if(self.debug):
			self.timesolve-=clock.time()
			casalog.post("DEBUG: start solve")
			
		self.getGains(solver,accumData,accumModel,accumWeight,accumFlag,goodbl)
				
		if(self.debug):
			self.timesolve+=clock.time()
			casalog.post("DEBUG: end solve")
			casalog.post("DEBUG: Time on solver- %.3f"%solver.timeOnSolver)
			casalog.post("DEBUG: Time on flagging- %.3f"%solver.timeOnFlaging)
					
		self.nflaginit+=solver.nflagsinit
		self.gainsOld=np.copy(self.gains)					
		self.snr=np.zeros(self.gains.shape)
		self.snr[self.antFlags]=np.abs(self.gains[self.antFlags])/np.abs(self.gains_er[self.antFlags])
		totalSol=self.gains.shape[0]*self.gains.shape[1]*self.gains.shape[2]
		self.nflagfinal+=np.sum(self.antFlags)
		self.ntotalsol+=totalSol

		casalog.post("%d/%d : %d out of %d solutions flagged"%(self.isol+1,self.solintMap[ispw].shape[0],totalSol-np.sum(self.antFlags),totalSol))
			
	def accumulateAndSolvePerSPW(self):		
		self.gainsOld=np.copy(self.gains)
		ispw=self.ispw
		dosolve=True
		accumScan=-1
		while(dosolve):
			
			dosolve=self.accum[ispw].nextIter()
			if(dosolve):					
				accumTime,accumScan,accumSPW,accumField,accumData,accumModel,accumWeight,accumFlag=self.accum[ispw].getAccum()
				self.coreSolve(accumData,accumModel,accumWeight,accumFlag)
				
				self.putInTable(self.gains,self.gains_er,self.snr,np.logical_not(self.antFlags),accumTime/self.solintMap[ispw][self.isol],accumScan,accumField,accumSPW)
				self.isol+=1
				if(self.isol<len(self.solintMap[ispw])):
					self.accum[ispw].setnsolint(self.solintMap[ispw][self.isol])			
				self.accum[ispw].resetAccum()
		endflag=self.accum[ispw].getEndflag()
		#scan=[-1,-1]
		#if(endflag):		
		#	dataAll=self.ms.getdata(['SCAN_NUMBER'],ifraxis=True)	
		#	scan=dataAll['scan_number']

		return endflag,accumScan
		
	def accumulateAndSolveCombSPW(self):		
		self.gainsOld=np.copy(self.gains)
		
		dosolve=True
		accumScan=-1
		while(dosolve):
			for ispw in range(len(self.spwIndex)):
				dosolve=self.accum[ispw].nextIter()
			
			if(dosolve):					
				accumTime,accumScan,accumSPW,accumField,accumData,accumModel,accumWeight,accumFlag=self.accum[0].getAccum()
				for ispw in range(1,len(self.spwIndex)):
					accumTime,accumScan,accumSPW,accumField,accumData_,accumModel_,accumWeight_,accumFlag_=self.accum[ispw].getAccum()
					accumData=np.append(accumData,accumData_,axis=2)
					accumModel=np.append(accumModel,accumModel_,axis=2)
					accumWeight=np.append(accumWeight,accumWeight_,axis=2)
					accumFlag=np.append(accumFlag,accumFlag_,axis=2)
				
				self.coreSolve(accumData,accumModel,accumWeight,accumFlag)
				
				self.isol+=1
				for ispw in range(len(self.spwIndex)):
					self.putInTable(self.gains,self.gains_er,self.snr,np.logical_not(self.antFlags),accumTime/self.solintMap[ispw][self.isol-1],accumScan,accumField,self.spwIndex[ispw])
					if(self.isol<len(self.solintMap[ispw])):
						self.accum[ispw].setnsolint(self.solintMap[ispw][self.isol])			
					self.accum[ispw].resetAccum()
		endflag=self.accum[0].getEndflag()

		#scan=[-1,-1]
		#if(endflag):		
		#	dataAll=self.ms.getdata(['SCAN_NUMBER'],ifraxis=True)	
		#	scan=dataAll['scan_number']
		
		return endflag,accumScan


	def findDataShape(self):
		self.attachMS()		
		self.ms.iterinit(interval=self.tblock,adddefaultsortcolumns=False,maxrows=1) 
                self.ms.iterorigin()                         
                dataAll=self.ms.getdata([self.dcolumn.upper()],ifraxis=True) 
		d=dataAll[self.dcolumn.lower()]
		self.nChan=d.shape[1]
		self.nBaseline=d.shape[2]
		self.ms.iterend()
 
	def solve(self):			
		self.hasmodel=False
		self.hasmodelinit=False
		    
		self.keystoget=[self.dcolumn.upper(),'FLAG','FIELD_ID','DATA_DESC_ID',
						'ANTENNA1','ANTENNA2','UVW','TIME','SCAN_NUMBER','SIGMA']
		
		#Initializing measurement set
		self.attachMS()
		
		#Doing initial iteration to find nCorr,ant1,ant2
		self.ms.iterinit(interval=self.tblock,adddefaultsortcolumns=False,maxrows=1) 
                self.ms.iterorigin()  
                if(self.debug):
                	casalog.post("DEBUG: Before getdata")                               
                dataAll=self.ms.getdata([self.dcolumn.upper(),'SCAN_NUMBER','MODEL_DATA'],ifraxis=True) 
		if(self.debug):
               		casalog.post("DEBUG: Getting dimensions")      
                if(len(dataAll['model_data'])!=0):
                	self.hasmodel=True
                	self.hasmodelinit=True
                	self.keystoget.append('MODEL_DATA')      	
           
		d=dataAll[self.dcolumn.lower()]
		
		scanstop_=np.min(dataAll['scan_number'])-1 #The -1 to ensure that calibration on the first block is done properly
		scanstop=np.max(dataAll['scan_number'])
		#Initializing accumulators and gains
		if(self.debug):
			casalog.post("DEBUG: Initializing accumulators")
		

		
		#Function to solve on blockID 
		self.ms.iterend()
		self.ms.close()
		self.attachMS()
		self.ms.iterinit(interval=self.tblock,adddefaultsortcolumns=False) 				
			
				
		blocki=0	
		
		scanrange=''
		if(len(self.gaintable)!=0):
			if(scanstop_+1==scanstop):
				scanrange='%d'%scanstop
			else:
				scanrange='%d~%d'%(scanstop_+1,scanstop)
			self.corrupt(scanrange)
			#Reattaching measurement set to ensure model_data column is detected
			if(not self.hasmodelinit):
				self.ms.close()
				self.attachMS()	
		
			
		if(not self.combinespw):
			for ispw in range(len(self.spwIndex)):

				self.accum[ispw].setHasModel(self.hasmodel)
				self.ispw=ispw
				self.isol=0
				self.findDataShape()
				self.initializeGains()
				
				self.accum[self.ispw].setnsolint(self.solintMap[self.ispw][0])
				endflag=True
				while(endflag):	
					scanstop_=scanstop
					endflag,scanstop=self.accumulateAndSolvePerSPW()	
					if(scanstop_<scanstop):
						self.restoremodel(scanrange)
						if(scanstop_+1==scanstop):
							scanrange='%d'%scanstop
						else:
							scanrange='%d~%d'%(scanstop_+1,scanstop)
						self.corrupt(scanrange)	
			
		else:	
			self.initializeGains()
			self.findDataShape()
			for ispw in range(len(self.spwIndex)):
				self.accum[ispw].setHasModel(self.hasmodel)
				#self.ispw=ispw

				self.accum[ispw].setnsolint(self.solintMap[ispw][0])

			self.isol=0
			endflag=True
			while(endflag):	
				scanstop_=scanstop
				endflag,scanstop=self.accumulateAndSolveCombSPW()			
				if(scanstop_<scanstop):
					self.restoremodel(scanrange)
					if(scanstop_+1==scanstop):
						scanrange='%d'%scanstop
					else:
						scanrange='%d~%d'%(scanstop_+1,scanstop)
					self.corrupt(scanrange)
	
			
		
		self.restoremodel(scanrange)
		self.ms.iterend()
		self.ms.close()			
		if(not self.hasmodelinit and self.hasmodel):
			delmod.delmod(vis=self.vis,scr=True)
		if(self.debug):
			casalog.post("DEBUG: Total time in MS I/O: %.3f"%self.timemsio)
			casalog.post("DEBUG: Total time in accumulation: %.3f"%self.timeaccum)
			casalog.post("DEBUG: Total time in mean/median computation: %.3f"%self.timemean)
			casalog.post("DEBUG: Total time in solve: %.3f"%self.timesolve)
		if(self.combinepoln):
			self.nflaginit*=2
		self.nflaginit=self.ntotalsol-self.nflaginit
		self.nflagfinal=self.ntotalsol-self.nflagfinal
		nflagrobust=self.nflagfinal-self.nflaginit
		casalog.post("Solution statistics:")
		casalog.post("%d/%d (%.2f %%) solutions were already flagged in ms"%(self.nflaginit,self.ntotalsol,(100.0*self.nflaginit)/self.ntotalsol))
		casalog.post("%d/%d (%.2f %%) solutions flagged due to low SNR and/or baseline count"%(nflagrobust,self.ntotalsol,(100.0*nflagrobust)/self.ntotalsol))
		casalog.post("%d/%d (%.2f %%) solutions flagged in total"%(self.nflagfinal,self.ntotalsol,(100.0*self.nflagfinal)/self.ntotalsol))		

		return 1
