import calRcore.coreIO
import numpy as np
import re
from taskinit import casalog
class bandpassR(calRcore.coreIO.calibSolver):
	def __init__(self,dividebychanzero=False,chanzerorange='',normamp=True,
		     zerophase=True,preaverage=False,
                     solnorm=True,normchanrange='',normampgains=True,
                     zerophasegains=True,**kwargs):
                     
		calRcore.coreIO.calibSolver.__init__(self,**kwargs)
		self.createCaltable("B Jones",singleChan=False)
		#self.correctCaltableSPW()
		self.dividebychanzero=dividebychanzero
		self.normamp=normamp
		self.zerophase=zerophase
		self.preaverage=preaverage
		self.chanzerorangetxt=chanzerorange
		
		self.solnorm=solnorm
		self.normchanrangetxt=normchanrange
		self.normampgains=normampgains
                self.zerophasegains=zerophasegains
		
		self.combinepoln=False
		self.combinespw=False
	'''	
	def correctCaltableSPW(self):
		if(len(self.chanindex)==0):
			return
		tb=self.tb
		slce=np.s_[self.chanindex[0,1]:self.chanindex[0,2]:self.chanindex[0,3]]
		tb.open(self.caltable+"/SPECTRAL_WINDOW",nomodify=False)
		chan_freq=tb.getcol('CHAN_FREQ')
		
		chan_freq=chan_freq[slce]
		tb.putcol('CHAN_FREQ',chan_freq)
		
		
		chan_width=tb.getcol('CHAN_WIDTH')
		chan_width=chan_width[slce]
		tb.putcol('CHAN_WIDTH',chan_width)
		
		effective_bw=tb.getcol('EFFECTIVE_BW')
		effective_bw=effective_bw[slce]
		tb.putcol('EFFECTIVE_BW',effective_bw)
		
		resolution=tb.getcol('RESOLUTION')
		resolution=resolution[slce]
		tb.putcol('RESOLUTION',resolution)
		
		tb.putcol('NUM_CHAN',np.array([len(chan_freq)]))
		
		tb.putcol('TOTAL_BANDWIDTH',np.array([np.sum(chan_width)]))
		
		tb.close()
	'''
	def parseRange(self,chanzerorange):
		matchObj=re.match("(.*)~(.*)",chanzerorange)			
		if(matchObj!=None and len(matchObj.groups())==2):
			try:
				rnge=[int(matchObj.group(1)),int(matchObj.group(2))+1]
				if(rnge[0]>rnge[1]):
					error=True
			except:
				error=True		
				
		else:
			matchObj=re.match("<(.*)",chanzerorange)
			if(matchObj!=None and len(matchObj.groups())==1):
				try:
					rnge=[0,int(matchObj.group(1))]
				except:
					error=True			
			else:
				matchObj=re.match("(.*)>",chanzerorange)
				if(matchObj!=None and len(matchObj.groups())==1):
					try:
						rnge=[int(matchObj.group(1)),self.nChan]
					except:
						error=True			
				else:
					matchObj=re.match("(.*)",chanzerorange)
					if(matchObj!=None and len(matchObj.groups())==1):
						try:
							rnge=[int(matchObj.group(1)),int(matchObj.group(1))+1]
						except:
							error=True
					else:
						error=True

		return rnge
		
	def initializeGains(self):		
		self.gains=np.zeros((self.nCorr,self.nChan,self.Nant),dtype=np.complex128) 
		self.gains_er=np.zeros((self.nCorr,self.nChan,self.Nant))
		self.antFlags=np.zeros((self.nCorr,self.nChan,self.Nant),dtype=np.bool)

		if(self.dividebychanzero):
			if(self.chanzerorangetxt==''):
				self.chanzerorange=np.array([0.25*self.nChan,0.75*self.nChan],dtype=np.int)
			else:
				self.chanzerorange=self.parseRange(self.chanzerorangetxt)
				
		if(self.solnorm):
			if(self.normchanrangetxt==''):
				self.normchanrange=np.array([0.25*self.nChan,0.75*self.nChan],dtype=np.int)
			else:
				self.normchanrange=self.parseRange(self.normchanrangetxt)
#				if(len(self.chanindex)!=0):
#					self.normchanrange-=self.chanindex[0,1]
			
			

	
	def dividechanzero(self,accumd,accumfl):
		samples=self.solintMap[self.soli]
		

		if(self.debug):
			casalog.post("DEBUG: chanzerorange: %d,%d"%(self.chanzerorange[0],self.chanzerorange[1]))	
		chanzero=accumd[:samples,:,self.chanzerorange[0]:self.chanzerorange[1],:]
		chanzeroflags=accumfl[:samples,:,self.chanzerorange[0]:self.chanzerorange[1],:]
		chanzero=self.centralValue(accumd[:samples,:,self.chanzerorange[0]:self.chanzerorange[1],:],flags=chanzeroflags,axis=1,keepdims=True)
			
		chanzeroflags=np.mean(chanzeroflags,axis=1,keepdims=True)
		chanzeroflags=chanzeroflags>0.1
			
		if(self.preaverage):
			chanzero=self.centralValue(chanzero,flags=chanzeroflags,weights=None,axis=3,keepdims=True)	

			
		if(self.normamp and self.zerophase):
			accumd/=chanzero
		elif(self.normamp):
			accumd/=np.abs(chanzero)
		elif(self.normphase):
			chanzero/=np.abs(chanzero)
			accumd/=chanzero
		return accumd
		'''
		if(self.debug):	
			print("DEBUG: computing central value")
		d=self.centralValue(self.accumd[:,:,:,:samples],flags=self.accumfl[:,:,:,:samples],weights=None,axis=3)
		
		wt=np.sqrt(np.mean(self.accumfl[:,:,:,:samples],axis=3))
		fl_tmp=np.mean(self.accumfl[:,:,:,:samples],axis=3)
		fl=fl_tmp>0.1		
		if(self.debug):
			print("DEBUG: done computing central value")
		return d,wt,fl
		'''
	
	def normalizeGains(self):
		if(self.debug):
				print("DEBUG: normchanrange : %d,%d"%(self.normchanrange[0],self.normchanrange[1]))	
		chanzero=np.mean(self.gains[:,self.normchanrange[0]:self.normchanrange[1],:],axis=1,keepdims=True)
		if(self.normampgains and self.zerophasegains):
			self.gains/=chanzero
			self.gains_er/=np.abs(chanzero)
		elif(self.normampgains):
			self.gains/=np.abs(chanzero)
			self.gains_er/=np.abs(chanzero)
		elif(self.zerophasegains):
			chanzero/=np.abs(chanzero)
			self.gains/=chanzero
		
			
	def getGains(self,solver,accumd,accummodel,accumwt,accumfl,goodbl):
		if(accumd.shape[2]!=self.gains.shape[1]):
			print("Sub-selection of channels not supported in bandpassR")
			self.error=True
			return
		if(self.dividebychanzero):
			accumd=self.dividechanzero(accumd,accumfl) #for divide by chanzero	
		nsamples=self.solintMap[self.ispw][self.isol]
		for ichan in range(0,self.nChan):
			if(self.debug):
				casalog.post('DEBUG: solving channel %d'%ichan)
				
			for icorr in range(0,self.nCorr):
				self.gains[icorr,ichan,:],self.gains_er[icorr,ichan,:],self.antFlags[icorr,ichan,:]=solver.solve(accumd[:nsamples,goodbl,ichan:ichan+1,icorr],accummodel[:nsamples,goodbl,ichan:ichan+1,icorr],accumfl[:nsamples,goodbl,ichan:ichan+1,icorr],accumwt[:nsamples,goodbl,ichan:ichan+1,icorr])
				if(not self.antFlags[icorr,ichan,self.refant]):
					casalog.post("change in refant")
					self.gains[icorr,ichan,:]=self.gains[icorr,ichan,:]*np.exp(1j*np.angle(self.gainsOld[icorr,ichan,solver.refant]))					
					solver.refant=self.refant
		if(self.solnorm):
			self.normalizeGains()
		if(self.debug):
			print ''
