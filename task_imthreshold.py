
import remove_negetive



def imthreshold(imname=None,boxsize=1,threshold=0,mtmfs=True,nterm=3,debugmode=False):

	remove_negetive.removenegetive(imname=imname,boxsize=boxsize,threshold=threshold,mtmfs=mtmfs, nterms=nterm,plotslices=False,debugmode=debugmode)
		
	return None

