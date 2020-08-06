using namespace std;
#include <complex>
#include <iostream>
#include <memory.h>
extern "C" {
void updateGains(complex<double> *X,complex<double> *g,int *ant1,int *ant2, double *wt,
                 bool *flags, int Nb,int Nant,double alpha,int refant)
{

    complex<double> *num=new complex<double>[Nant];
    complex<double> *denom=new complex<double>[Nant];
    complex<double> *update=new complex<double>[Nant];
    double* gsqr=new double[Nant];
    complex<double> ref;
    int i,j;
    for(i=0; i<Nant;i++) {
        num[i]=0.0;
        denom[i]=0.0;
        gsqr[i]=pow(abs(g[i]),2.00);

    }
  
    for(int k=0; k<Nb;k++) {
    	if(flags[k]==0)
    		continue;
    	i=ant1[k];
    	j=ant2[k];
    	num[i]+=X[k]*g[j]*wt[k];
    	num[j]+=conj(X[k])*g[i]*wt[k];
    	denom[i]+=gsqr[j]*wt[k];
    	denom[j]+=gsqr[i]*wt[k];

    }
   
    //updating gains
    for(int i=0; i<Nant;i++) {
    	if(denom[i]==0.0)
    		update[i]=0.0;
    	else
    	        update[i]=num[i]/denom[i];
	g[i]=g[i]*(1.0-alpha)+alpha*(update[i]);
	
    }
    //setting relative to refant
    ref=polar(1.0,-arg(g[refant]));
    for(int i=0; i<Nant;i++) {
    	if(denom[i]!=0.0) 
		g[i]=g[i]*ref;	
    }

    delete[] gsqr;
    delete[] num;
    delete[] denom;
    delete[] update;

}

void updatePhase(complex<double> *X,complex<double> *g,int *ant1,int *ant2, double *wt,
                 bool *flags, int Nb,int Nant,double alpha,int refant)
{
    
    complex<double> *update=new complex<double>[Nant];
   //set update to zero!!
    complex<double> ref;
    int i,j;
  
   for(int k=0; k<Nb;k++) {

    	if(flags[k]==0)
    		continue;
    	i=ant1[k];
    	j=ant2[k];
    	update[i]+=X[k]*wt[k]/conj(g[j]);
    	update[j]+=conj(X[k])*wt[k]/conj(g[i]);
    }
    //updating gains
    for(int i=0; i<Nant;i++) {
    	
	g[i]=polar(1.0,arg(update[i]));
	
    }
    //setting relative to refant
    ref=polar(1.0,-arg(g[refant]));
    for(int i=0; i<Nant;i++) {
	g[i]=g[i]*ref;	
    }

    delete[] update;

}

int checkConv(complex<double> *gold,complex<double> *g,bool *antflags,int Nant,double tol)
{
	int cnt=0;
	for(int i=0;i<Nant;i++)
	{
		if(!antflags[i])
			continue;
		if(abs(real(g[i]-gold[i]))>tol*abs(g[i]))
			cnt++;
		if(abs(imag(g[i]-gold[i]))>tol*abs(g[i]))
			cnt++;
	}
	return cnt;
}

void getGains(complex<double> *X,complex<double> *g,int *ant1,
            int *ant2,double *wt,bool *antflags,bool *flags, int Nb,int Nant
            ,double alpha,int refant,double tol,bool doPhase,bool hasConverged)
{
    complex<double>* gold=new complex<double>[Nant];
    int nitermax=1000;
    //alpha=0.5;
    if(doPhase)
    	     updatePhase(X,g,ant1,ant2,wt,flags,Nb,Nant,alpha,refant);
    else
    	     updateGains(X,g,ant1,ant2,wt,flags,Nb,Nant,alpha,refant);  
    int niter=0;
    while(checkConv(gold,g,antflags,Nant,tol) && niter<nitermax)
    {
    	/**
	if(niter==2)
		alpha=0.75;
        if(niter>2)
        	alpha=0.9;
        if(niter>10)
        	alpha=0.5;
        **/
    	memcpy(gold,g,sizeof(g[0])*Nant);

	    if(doPhase)
	    	     updatePhase(X,g,ant1,ant2,wt,flags,Nb,Nant,alpha,refant);
	    else
	    	     updateGains(X,g,ant1,ant2,wt,flags,Nb,Nant,alpha,refant);  
	    	     
    	niter++;	
    	
    }
	if(niter>=nitermax)
    		hasConverged=0;
	else
		hasConverged=1;
    delete[] gold;
}
}
