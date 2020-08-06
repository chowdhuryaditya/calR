
import gaincalRIO



def gaincalR(vis=None,caltable=None,field='',spw='',selectdata=True,uvrange='',scan='',observation='',
	             solint='inf',combine='',refant=0,minblperant=6,
                     minsnr=5,calmode='ap',
                     robust=True,outlieralgorithm='RMS',
                     threshold=[7,5,4,3],dosolintflag=False,mingain=1e-2,combinepoln=False,combinespw=False,
                     solnorm=False,combinecorrnorm=True,usemedian=True,debugmode=False):

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
		dosolintflag=False
		threshold=[1e20]
		mingain=0.0
		
	c=gaincalRIO.gaincalR(vis=vis,caltable=caltable,field=field,spw=spw,uvrange=uvrange,scan=scan,
			observation=observation,solint=solint,combine=combine,refant=refant,minblperant=minblperant,
			minsnr=minsnr,calmode=calmode,outlieralgorithm=outlieralgorithm,
			threshold=threshold,dosolintflag=dosolintflag,mingain=mingain,
			gaintable=gaintable,gainfield=gainfield,interp=interp,spwmap=spwmap,
			combinepoln=combinepoln,combinespw=combinespw,
			solnorm=solnorm,combinecorrnorm=combinecorrnorm,usemedian=usemedian,
			debugmode=debugmode)
	
	c.solve()
	if(solnorm):
		c.normalize()
	return None

