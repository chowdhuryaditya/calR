import numpy as np
import matplotlib.pyplot as plt
from __casac__ import *
from taskinit import casalog
import os
ia=image.image()
def removenegetive(imname,boxsize,threshold,mtmfs=True,nterms=3,plotslices=False,debugmode=False):
	if(mtmfs):
		ia.open(imname+".model.tt0")
		im=ia.getchunk()
		ia.close()
	else:
		ia.open(imname+".model")
		im=ia.getchunk()
		ia.close()

	niter=0
	casalog.post("Total flux before thresholding : %.1f mJy"%(np.sum(im)*1e3))
	if(debugmode):
		casalog.post("Total number of points below threshold:%d"%np.sum(im<threshold))
	while(np.min(im)<threshold):
		minpoint=np.unravel_index(np.argmin(im, axis=None), im.shape)
		minimumflux=im[minpoint]
		localslice=im[minpoint[0]-boxsize//2:minpoint[0]+boxsize//2,minpoint[1]-boxsize//2:minpoint[1]+boxsize//2,0,0]
		if(debugmode):
			casalog.post("Iteration - %d; minimum at (%d,%d): %.3f mJy"%(niter+1,minpoint[0],minpoint[1],minimumflux*1e3))
		if(plotslices):
			plt.figure()
			plt.imshow(localslice)
			plt.colorbar()
			plt.show()

		localslice[localslice<np.abs(minimumflux)]=0.0
		im[minpoint[0]-boxsize//2:minpoint[0]+boxsize//2,minpoint[1]-boxsize//2:minpoint[1]+boxsize//2,0,0]=localslice
		niter+=1
	casalog.post("Total flux after thresholding : %.1f mJy"%(np.sum(im)*1e3))
	zeromask=im<1e-9
	if(mtmfs):
		for i in range(0,nterms):
			os.system("cp -r %s.model.tt%d %s_noneg.model.tt%d"%(imname,i,imname,i))
			ia.open("%s.model.tt%d"%(imname,i))
			imterm=ia.getchunk()
			ia.close()	
			imterm[zeromask]=0.0
			ia.open("%s_noneg.model.tt%d"%(imname,i))
			ia.putchunk(imterm)
			ia.close()
	else:
		os.system("cp -r %s.model %s_noneg.model"%(imname,imname))
		ia.open("%s.model"%(imname))
		imterm=ia.getchunk()
		ia.close()	
		imterm[zeromask]=0.0
		ia.open("%s_noneg.model"%(imname))
		ia.putchunk(imterm)
		ia.close()
		
#removenegetive("images/selfcal_round4_nterms4_briggs-0.5_cell0.6_imsize9k",150,0.0,plotslices=False)
#removenegetive("images/selfcal_round1_nterms3_briggs0.5_cell0.95_imsize6k",150,0.0,plotslices=False)
#removenegetive("peelingTests/images_nobs/round2",200,0.0,plotslices=False)
