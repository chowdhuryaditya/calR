import calRcore.coreIO
import numpy as np
from taskinit import casalog
class gaincalR(calRcore.coreIO.calibSolver):
	def __init__(self,combinepoln=False,combinespw=False,solnorm=False,combinecorrnorm=True,usemedian=True,**kwargs):
		calRcore.coreIO.calibSolver.__init__(self,**kwargs)
		self.createCaltable("G Jones",singleChan=True)
		self.solnorm=solnorm
		self.usemedian=usemedian
		self.combinepoln=combinepoln
		self.combinespw=combinespw
		self.combinecorrnorm=combinecorrnorm

		#if(self.combinespw):
		#	self.error=self.validateSolint()
		#	if(not self.error):
		#		print("Cannot combine Spectral Windows. Were different windows observed simultaneously?")
				
	def initializeGains(self):	
		if(self.nCorr==1):	
			self.gains=np.zeros((2,1,self.Nant),dtype=np.complex128) 
			self.gains_er=np.zeros((2,1,self.Nant))
			self.antFlags=np.zeros((2,1,self.Nant),dtype=np.bool)
		else:
			self.gains=np.zeros((self.nCorr,1,self.Nant),dtype=np.complex128) 
			self.gains_er=np.zeros((self.nCorr,1,self.Nant))
			self.antFlags=np.zeros((self.nCorr,1,self.Nant),dtype=np.bool)
	
	def validateSolint(self):
		nsol=np.sum(self.solintMap,axis=1)
		if(np.unique(nsol).shape[0]>1):
			return False
		if(np.sum(np.mean(self.solintMap,axis=0)!=self.solintMap[0])>0):
			return False
		return True  	
			
	def normalize(self):
		self.tb.open(self.caltable,nomodify=False)
		cparam=self.tb.getcol('CPARAM')
		paramerr=self.tb.getcol('PARAMERR')
		flags=self.tb.getcol('FLAG')
		if(self.usemedian):
			meanfunc=np.median
		else:
			meanfunc=np.mean
		if(self.combinecorrnorm):
			globalVal=meanfunc(np.abs(cparam[np.logical_not(flags)]))
			casalog.post("Central value of gains: %.3f"%globalVal)
			cparam/=globalVal
			paramerr/=globalVal
		else:
			for icorr in range(0,self.nCorr):
				globalVal=meanfunc(np.abs(cparam[icorr,:,:][np.logical_not(flags[icorr,:,:])]))
				cparam[icorr,:,:]/=globalVal
				paramerr[icorr,:,:]/=globalVal
		self.tb.putcol('CPARAM',cparam)
		self.tb.putcol('PARAMERR',paramerr)
		self.tb.putcol('FLAG',flags)
		self.tb.close()
				
	def getGains(self,solver,accumd,accummodel,accumwt,accumfl,goodbl):	
		nsamples=self.solintMap[self.ispw][self.isol]
		if(not self.combinepoln):	
		
			for j in range(0,self.nCorr):
				thisflags=accumfl[:nsamples,goodbl,:,j]
				if(np.sum(thisflags)==0):
					casalog.post("No unflagged samples")
					self.gains[j,0,:]=1+1j*0
					self.antFlags[j,0,:]=np.zeros(self.Nant,dtype=np.bool)
					continue
				self.gains[j,0,:],self.gains_er[j,0,:],self.antFlags[j,0,:]=solver.solve(accumd[:nsamples,goodbl,:,j],accummodel[:nsamples,goodbl,:,j],accumfl[:nsamples,goodbl,:,j],accumwt[:nsamples,goodbl,:,j])
				if(not self.antFlags[j,0,self.refant]):
					self.gains[j,0,:]=self.gains[j,0,:]*np.exp(1j*np.angle(self.gainsOld[j,0,solver.refant]))					
					solver.refant=self.refant
					
		else:
			thisflags=accumfl[:nsamples,goodbl,:,:]
			if(np.sum(thisflags)==0):
				casalog.post("No unflagged samples")
				self.gains[:,0,:]=1+1j*0
				self.antFlags[:,0,:]=np.zeros(self.Nant,dtype=np.bool)
				return

			nCorr=accumd.shape[3]
			nChan=accumd.shape[2]
			nBaseline=accumd.shape[1]
			newshape=(nsamples,nBaseline,nCorr*nChan)
			thisdata=np.swapaxes(accumd[:nsamples,:,:,:].reshape(newshape),0,1)
			thismodel=np.swapaxes(accummodel[:nsamples,:,:,:].reshape(newshape),0,1)
			thisflag=np.swapaxes(accumfl[:nsamples,:,:,:].reshape(newshape),0,1)
			accumwt_tile=np.tile(accumwt,(1,1,nChan//accumwt.shape[2],1))
			thiswt=np.swapaxes(accumwt_tile[:nsamples,:,:,:].reshape(newshape),0,1)
			self.gains[0,0,:],self.gains_er[0,0,:],self.antFlags[0,0,:]= solver.solve(thisdata[goodbl,:,:],thismodel[goodbl,:,:],thisflag[goodbl,:,:],thiswt[goodbl,:,:])
			if(not self.antFlags[0,0,self.refant]):
				self.gains[0,0,:]=self.gains[0,0,:]*np.exp(1j*np.angle(self.gainsOld[0,0,solver.refant]))					
				solver.refant=self.refant
			self.gains[1,:,:]=self.gains[0,:,:]
			self.gains_er[1,:,:]=self.gains_er[0,:,:]
			self.antFlags[1,:,:]=self.antFlags[0,:,:]
		
