import bandpassRIO

def bandpassR(vis=None,caltable=None,field='',spw='',selectdata=True,uvrange='',scan='',observation='',
	             solint='inf',combine='',refant=0,minblperant=6,
                     minsnr=5,calmode='ap',robust=True,outlieralgorithm='RMS',
                     threshold=[7,5,4,3],dosolintflag=True,mingain=1e-2,
                     dividebychanzero=False,chanzerorange='',normamp=True,zerophase=True,preaverage=False,
                     solnorm=True,normchanrange='',normampgains=True,zerophasegains=True,debugmode=False):

        if(not selectdata):
        	uvrange=''
        	scan=''
        	observation=''
        caltables=False
        if(not caltables):
        	gaintable=[]
        	gainfield=[]
        	interp=[]
        	spwmap=[]

	if(not robust):
		dosolintmedian=False
		threshold=[1e20]
		mingain=0.0
		
	if(not dividebychanzero):
		chanzerorange=''
		normamp=False
		zerophase=False
	if(not solnorm):
		normchanrange=''
		normampgains=False
		zerophasegains=False
	c=bandpassRIO.bandpassR(vis=vis,caltable=caltable,field=field,spw=spw,uvrange=uvrange,scan=scan,
			observation=observation,solint=solint,combine=combine,refant=refant,minblperant=minblperant,
			minsnr=minsnr,calmode=calmode,outlieralgorithm=outlieralgorithm,
			threshold=threshold,dosolintflag=dosolintflag,mingain=mingain,
			gaintable=gaintable,gainfield=gainfield,interp=interp,spwmap=spwmap,
			dividebychanzero=dividebychanzero,chanzerorange=chanzerorange,
			normamp=normamp,zerophase=zerophase,preaverage=preaverage,
                        solnorm=solnorm,normchanrange=normchanrange,normampgains=normampgains,zerophasegains=zerophasegains,debugmode=debugmode)
	if(c.error==True):
		return None
	
	c.solve()
	return None

