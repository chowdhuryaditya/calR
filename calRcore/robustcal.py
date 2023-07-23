import numpy as np
import matplotlib.pyplot as plt
import ctypes
import dummy
from ctypes.util import find_library
import os
import time as clock
from taskinit import casalog
pathname=os.path.dirname(dummy.__file__)
solver = ctypes.CDLL(pathname+"/calRcore/solver.so")

class rantsol:	
	def __init__(self,threshold,mingain,refant,minblperant,Nant,ant1,ant2,outlieralgorithm,dorobustmean,calmode,minsnr,dosolintflag,debug=False,delta=1e-7,alpha=0.3):
		self.delta=delta
		self.Nant=Nant
		self.Nb=ant1.shape[0]
		self.alpha=alpha
		self.ant1=ant1
		self.ant2=ant2
		self.refant=refant
		self.minblperant=minblperant
		self.nsigmaList=threshold
		self.outlieralgorithm=outlieralgorithm
		self.dorobustmean=dorobustmean
		self.dosolintflag=dosolintflag
		self.calmode=calmode
		self.mingain=mingain
		self.timeOnSolver=0.0
		self.timeOnFlaging=0.0
		self.nflagsinit=0
		self.minsnr=minsnr
		self.counterglobal=0
		self.counterlocal=0
		self.debug=debug
		self.niterFlagMax=100
		if(not self.debug):
			np.seterr(divide='ignore', invalid='ignore')
			np.warnings.filterwarnings('ignore')
	def outlierDetect(self,x,flags=[]):
		real=np.real(x)
		imag=np.imag(x)
		if(len(flags)>0):
			if(np.sum(flags)==0):
				return flags
			realflags=real[flags]
			imagflags=imag[flags]
		else:
			realflags=real
			imagflags=imag
		'''
		print(np.sum(flags))
		plt.figure()
		plt.plot(real[flags])
		plt.plot(imag[flags])
		plt.savefig("series_%d_%d.png"%(self.counterglobal,self.counterlocal))
		self.counterlocal+=1
		'''
		
		if(self.outlieralgorithm=='RMS'):
			sigma_r=np.sqrt(np.mean(realflags**2))
			sigma_i=np.sqrt(np.mean(imagflags**2))
		else:
			sigma_r=1.4826*np.median(np.abs(realflags))	
			sigma_i=1.4826*np.median(np.abs(imagflags))	

		flags_r=np.abs(real)<=self.nsigma*sigma_r	
		flags_i=np.abs(imag)<=self.nsigma*sigma_i
		return np.logical_and(flags_r,flags_i)

	def getFlags(self,X,flags,olderflags,antflags,g):
		self.timeOnFlaging-=clock.time()
		gs=np.conj(g)		
		fl=np.copy(flags)
		antlist=np.arange(self.Nant)
		if(self.calmode=='p'):
			fl[flags]=self.outlierDetect(X[flags]/np.abs(X[flags])-g[self.ant1][flags]*gs[self.ant2][flags],olderflags[flags])
		else:
			fl[flags]=self.outlierDetect(X[flags]-g[self.ant1][flags]*gs[self.ant2][flags],olderflags[flags])
		antflags[np.abs(g)<self.mingain]=False
		for i in range(0,self.Nant):	
			if((np.sum(fl[self.ant1==i])+np.sum(fl[self.ant2==i]))<self.minblperant):
				antflags[i]=False
		for i in range(0,self.Nant):
			if(not antflags[i]):
				fl[self.ant1==i]=False
				fl[self.ant2==i]=False
		self.timeOnFlaging+=clock.time()
		return antflags,fl
	'''
	def getFlags2(self,X,flags,olderflags,antflags,g):
		self.timeOnFlaging-=clock.time()
		gs=np.conj(g)		
		fl=np.copy(flags)
		Xs=np.conj(X)
		antlist=np.arange(self.Nant)
		#if(self.calmode=='p'):
		#	fl[flags]=self.outlierDetect(X[flags]/np.abs(X[flags])-g[self.ant1][flags]*gs[self.ant2][flags],olderflags[flags])
		#else:
		#	fl[flags]=self.outlierDetect(X[flags]-g[self.ant1][flags]*gs[self.ant2][flags],olderflags[flags])
		antflags[np.abs(g)<self.mingain]=False
		for i in range(0,self.Nant):
			index1=(self.ant1==i)
			index2=(self.ant2==i)	
			index=np.logical_or(index1,index2)
			
			g_i=np.append(X[index1]/gs[self.ant2[index1]],Xs[index2]/gs[self.ant1[index2]])
			flags_i=fl[index]
			fl[index][flags_i]=self.outlierDetect((g_i-g[i])[flags_i],olderflags[index][flags_i])
		for i in range(0,self.Nant):
			#print("%d: %d"%(i,np.sum(fl[self.ant1==i])+np.sum(fl[self.ant2==i])))
			if((np.sum(fl[self.ant1==i])+np.sum(fl[self.ant2==i]))<self.minblperant):
				antflags[i]=False
		
		for i in range(0,self.Nant):
			if(not antflags[i]):
				fl[self.ant1==i]=False
				fl[self.ant2==i]=False
		self.timeOnFlaging+=clock.time()
		return antflags,fl
	'''		
	def getGainSNRAIPS(self,X,flags,g): 
		#This, I think, is identical to what AIPS does.
		g_s=np.conj(g)
		gsqr=np.real(g*np.conj(g))
		g_er=np.zeros(self.Nant)		
		antIndx=np.arange(0,self.Nant)
		X_m=g_s[self.ant1]*g[self.ant2]
		for i in range(0,self.Nant):
			Xi=np.append(X[self.ant1==i],np.conj(X[self.ant2==i]))
			Xi_m=np.append(X_m[self.ant1==i],np.conj(X_m[self.ant2==i]))
			flagsi=np.append(flags[self.ant1==i],flags[self.ant2==i])

			g_er[i]=np.sqrt(np.mean(np.angle(Xi[flagsi]*Xi_m[flagsi])**2))
		return 1.0/g_er
		
	def getGainSNR(self,X,flags,g): 
		#A better variant?
		g_s=np.conj(g)
		gsqr=np.real(g*np.conj(g))
		g_er=1e10*np.ones(self.Nant)		
		antIndx=np.arange(0,self.Nant)
		
		X_m=g_s[self.ant1]*g[self.ant2]
		for i in range(0,self.Nant):			
			flagsi=np.append(flags[self.ant1==i],flags[self.ant2==i])
			if(np.sum(flagsi)==0):
				continue
			gi=np.append(X[self.ant1==i]/g_s[self.ant2[self.ant1==i]],np.conj(X[self.ant2==i])/g_s[self.ant1[self.ant2==i]])
			gsqri=np.append(gsqr[self.ant2[self.ant1==i]],gsqr[self.ant1[self.ant2==i]])
			g_er[i]=np.sqrt(np.average(np.abs(gi[flagsi]-g[i])**2,weights=gsqri[flagsi])/len(gi))
			
		return g_er
		
		
	def centralValue(self,x,model,flags,weights=[]):
		if(len(weights)==0):
			wf=flags
		else:
			wf=weights*flags
		fl=np.mean(np.mean(flags,axis=1),axis=1)
		X=np.sum(np.sum(x*np.conj(model)*wf,axis=1),axis=1)/np.sum(np.sum(model*np.conj(model)*wf,axis=1),axis=1)
		retflags=(fl>0.1)	
		#X=np.sum(np.sum(x*flags,axis=1),axis=1)/np.sum(np.sum(flags,axis=1),axis=1)
		return X,retflags
		
	def getweights(self,model,flags,weights):	
		wt=np.sum(np.sum(weights*np.abs(model)*flags,axis=1),axis=1)
		return wt
		
	def flagSubInt(self,X,model,flags,g,coarseflags):
		self.timeOnFlaging-=clock.time()
		flags[np.logical_not(coarseflags),:,:]=False
		gs=np.conj(g)	
		Xm=np.expand_dims(np.expand_dims(g[self.ant1]*gs[self.ant2],axis=1),axis=1)*model
		if(self.calmode=='p'):
			residuals=X/np.abs(X)-Xm/np.abs(Xm)
		else:
			
			residuals=X-Xm
		flags[flags]=self.outlierDetect(residuals[flags])
		self.timeOnFlaging+=clock.time()
		return flags
		
	def solve(self,X_full,model_full,flags_full,weights_full):
		#disable solint flag if there is only one sub-sample per solint
		X_full=np.complex128(X_full)
		model_full=np.complex128(model_full)
		weights_full=weights_full.astype("double")
		if(X_full.shape[1]*X_full.shape[2]==1):
			self.dosolintflag=False
		#avoiding infinity and NaN's
		indx=(np.abs(model_full)<1e-20)		
		flags_full[indx]=False	
		indx=(np.abs(X_full)<1e-20)		
		flags_full[indx]=False	
		#weights_full=1.0/sigma_full**2
		#weights_full[sigma_full<1e-20]=0.0
		weights_full[np.isnan(weights_full)]=0.0
		weights_full[np.isinf(weights_full)]=0.0
		X,flags=self.centralValue(X_full,model_full,flags_full,weights=weights_full)
		wt=self.getweights(model_full,flags_full,weights_full)
		self.counterlocal=0
		doPhase=False
		hasGainConverged=False
		fl=np.ones(self.Nb,dtype=np.bool)
		antflags=np.ones(self.Nant,dtype=np.bool)
		antlist=np.arange(self.Nant)
				
		if(self.calmode=='ap'):
			g_old=np.ones(self.Nant)+1j*np.zeros(self.Nant)
		elif(self.calmode=='p'):
			g_old=np.ones(self.Nant)+1j*np.zeros(self.Nant)
			doPhase=True
			#X/=np.abs(X)
		else:
			g_old=np.ones(self.Nant)
			X=np.abs(X)


		g_er=[]
		
		def updateGains():		
			
			self.timeOnSolver-=clock.time()
			solver.getGains(ctypes.c_void_p(X.ctypes.data),
					ctypes.c_void_p(g.ctypes.data),
					self.ant1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
					self.ant2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
					wt.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
					antflags.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
					,fl.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
					ctypes.c_int(self.Nb),
					ctypes.c_int(self.Nant),
					ctypes.c_double(self.alpha),ctypes.c_int(self.refant),
					ctypes.c_double(self.delta),ctypes.c_bool(doPhase),ctypes.c_bool(hasGainConverged))
			self.timeOnSolver+=clock.time()
		def checkConv():
			cnt=solver.checkConv(ctypes.c_void_p(g_old.ctypes.data),ctypes.c_void_p(g.ctypes.data),
					 antflags.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
					 ctypes.c_int(self.Nant),ctypes.c_double(self.delta))
			return cnt
		self.nsigma=1e20
		g=np.copy(g_old)
		antflags,fl=self.getFlags(X,flags,fl,antflags,g) 
		#if(self.debug):
		#	print("Initial flags: %d"%np.sum(antflags))
		self.nflagsinit+=np.sum(antflags)

		def chisqr():
			gs=np.conj(g)
			chi2=np.sqrt(np.sum(np.abs(X[flags]-g[self.ant1][flags]*gs[self.ant2][flags])**2))
			return chi2
		#initchi2=chisqr()
		converged=False
		while(not converged):
			for i in range(0,len(self.nsigmaList)):
				#print(np.mean(flags_full))
				self.nsigma=self.nsigmaList[i]
				if(i!=0):
					antflags,fl=self.getFlags(X,flags,fl,antflags,g) 
					g_old=np.copy(g)
					if(np.sum(antflags)==0):
						break
				updateGains()
				if(self.dosolintflag):
					flags_full=self.flagSubInt(X_full,model_full,flags_full,g,fl)	
					X,flags=self.centralValue(X_full,model_full,flags_full,weights=weights_full)
					fl[np.logical_not(flags)]=False
				#print(self.nsigmaList[i])
				niter=0
				while(checkConv()>0 and niter<self.niterFlagMax):
					g_old=np.copy(g)
					updateGains()
					antflags,fl=self.getFlags(X,flags,fl,antflags,g) 
					#print(np.sum(fl))
					#print(np.sum(antflags))
					niter+=1
					
					if(not antflags[self.refant]):
						casalog.post("Refant flagged, trying to find another refant.")
						if(np.sum(antflags)==0):
							casalog.post("Could not find another unflagged antenna; all solutions flagged")
							antflags=np.zeros(self.Nant,dtype=np.bool)
							break
							
						else:
							self.refant=antlist[antflags][0]
			#print("hasGainConverged: %d"%hasGainConverged)
			#finalchi2=chisqr()
			#print("fractional chisqr change:%f"%((initchi2-finalchi2)/initchi2))
			if(self.calmode=='p'):				
				g_er=self.getGainSNR(X/np.abs(X),fl,g)
			else:
				g_er=self.getGainSNR(X,fl,g)
			antflags_bk=np.copy(antflags)
			antflags[np.abs(g)/g_er<self.minsnr]=False
			if(np.sum(antflags_bk!=antflags)>0):
				converged=False
			else:
				converged=True
			if(np.sum(antflags)==0):
				break

		if(np.sum(antflags)>0 and hasGainConverged==True):
			print("Warning: Gain solutions have not converged; Maximum iteration exceeded")	
		g[np.logical_not(antflags)]=1.0+1j*0.0
		self.counterglobal+=1
		#if(self.debug):
		#	print("Final flags: %d"%np.sum(antflags))
		return g,g_er,antflags
	
