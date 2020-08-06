
import remove_negetive



def imthreshold(imname=None,boxsize=-1,threshold=0,mtmfs=True,nterm=3,debugmode=False):
	if(boxsize==1):
		print("ERROR: Parameter boxsize can either be -1 (i.e. imsize) or a positive integer >1")
		return None
	remove_negetive.removenegetive(imname=imname,boxsize=boxsize,threshold=threshold,mtmfs=mtmfs, nterms=nterm,plotslices=False,debugmode=debugmode)
		
	return None

