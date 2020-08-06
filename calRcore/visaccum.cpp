#include<visaccum.h>

VisAccum::VisAccum(string msname,string fieldname, string spw,string uvrange,string scan,string observation,string poln,int ntimesample,double timeInverval,int dataDescIndx,int spwIndx)
{
	msSelect(msname,fieldname,spw,uvrange,scan,observation,poln,dataDescIndx,spwIndx);
	MeasurementSet msOriginal(msname);
	nant=ROMSAntennaColumns(msOriginal.antenna()).nrow();
	
	initIter(timeInverval);
	initAccum(getDim(ntimesample));	
	resetAccum();
	tstart=0;
	hasmodel=false;
	
}
void VisAccum::setHasModel(bool _hasmodel)
{
	hasmodel=_hasmodel;
}
IPosition VisAccum::getDim(int ntimesample)
{
	IPosition dim(4);
	dim[2]=(nant*(nant-1))/2;
	
	dim[3]=ntimesample;
	ArrayColumn<Complex> _data (msIter->table(),  MS::columnName(MS::DATA));
	Array<Complex> data=_data.getColumn(cslice);

	dim[0]=data.shape()[0];
	dim[1]=data.shape()[1];
	
	return dim;
}

void VisAccum::msSelect(string &msname,string &fieldname, string &spw,string &uvrange,string &scan,string &observation,string &poln,int dataDescIndx,int spwIndx)
{
	MeasurementSet msOriginal(msname);

	string timeExpr,antennaExpr,fieldExpr,spwExpr, uvDistExpr,taQLExpr,polnExpr, scanExpr,arrayExpr,stateExpr,obsExpr,feedExpr;
	fieldExpr=fieldname;
	spwExpr=spw;
	uvDistExpr=uvrange;
	scanExpr=scan;
	obsExpr=observation;
	polnExpr=poln;
	taQLExpr="ANTENNA1!=ANTENNA2";
	MSSelection select;
	
	select.setFieldExpr(fieldExpr);
	select.setSpwExpr(spwExpr);
	select.setUvDistExpr(uvDistExpr);
	select.setScanExpr(scanExpr);
	select.setArrayExpr(arrayExpr);
	select.setObservationExpr(obsExpr);
	select.setPolnExpr(polnExpr);
	//cout<<fieldExpr<<":"<<spwExpr<<":"<<scanExpr<<":"<<arrayExpr<<":"<<obsExpr<<endl;
	Vector< Vector<Slice> > chanslices;
	Vector< Vector<Slice> > corrslices;
	
	ostringstream thisSpw;
	thisSpw<<spwIndx; 
	
	mssSetData2(msOriginal,ms,"",timeExpr,antennaExpr,fieldExpr,thisSpw.str(), uvDistExpr,taQLExpr,polnExpr, scanExpr,arrayExpr,stateExpr,obsExpr,feedExpr);

	select.getChanSlices(chanslices,&msOriginal);
	select.getCorrSlices(corrslices,&msOriginal);
	cslice.resize(2);
	cslice[0]=corrslices[0];
	cslice[1]=chanslices[dataDescIndx];
	//cout<<cslice<<endl;
}
Vector<int> VisAccum::getBaselineIndx(const Vector<int>& ant1,const Vector<int>& ant2)
{
	Vector<int> bindxArr(ant1.shape());
	int bindxData=0;
	int bindx=0;
	//cout<<ant1.shape()<<endl;
	//cout<<nant<<endl;
	for(int i=0;i<ant1.shape()[0];i++)
	{
		bindxArr[i]=(ant2[i]-ant1[i]-1)+(ant1[i]*((nant-1)+(nant-ant1[i])))/2;
		//cout<<i<<","<<ant1[i]<<","<<ant2[i]<<","<<bindxArr[i]<<endl;
	}
	/**
	for(int ii=0;ii<nant;ii++)
	{
		for(int jj=ii+1;jj<nant;jj++)
		{
			if(ant1[bindxData]==ii && ant2[bindxData]==jj)
			{
				cout<<ii<<","<<jj<<endl;
				bindxArr[bindxData]=bindx;
				bindxData++;
			}
			bindx++;
		}
	}
	cout<<bindx<<endl;

	**/
	return bindxArr;
}
Vector<Slice> VisAccum::getTimeStrides(const Vector<double>& time)
{
	int itime=0,i=0;
	double curt;
	Vector<Slice> tslice(time.shape()[0]); //excessively long!
	int startindxt=0;
	curt=time[startindxt];
	while(i<time.shape()[0])
	{
		
		if(time[i]>curt)
		{
			tslice[itime]=Slice(startindxt,i-startindxt,1);
			curt=time[i];
			itime++;
			startindxt=i;
		}
		i++;
		
	}
	tslice[itime]=Slice(startindxt,i-startindxt,1);
	return tslice(Slice(0,itime+1,1));
}
void VisAccum::initAccum(IPosition dim)
{
	accumData=Array<Complex>(dim);
	accumModel=Array<Complex>(dim);
	accumFlag=Array<bool>(dim);
	ArrayColumn<float> _weight; 
	try
	{
		hasWtSpec=true;
	 	_weight=ArrayColumn<float>(msIter->table(), MS::columnName(MS::SIGMA_SPECTRUM));
		IPosition tmp=_weight.shape(0); //CHECK IMPLICATIONS OF THIS LINE
		cout<<"testing this:"<<tmp<<endl;
	}
	catch(...)
	{
		hasWtSpec=false;
		_weight=ArrayColumn<float>(msIter->table(), MS::columnName(MS::SIGMA));
		dim[1] =1;
	}
	accumWeight=Array<float>(dim);

	iterData=new ArrayIterator<Complex>(accumData,3);
	iterWeight=new ArrayIterator<float>(accumWeight,3);
	iterModel=new ArrayIterator<Complex>(accumModel,3);
	iterFlag=new ArrayIterator<bool>(accumFlag,3);
}
void VisAccum::resetAccum()
{
	accumFlag=true;
	accumTime=0.0;
	accumScan=-1;
	accumSPW=-1;
	accumModel=1;
	accumWeight=(float)0.0;
	iterData->reset();
	iterWeight->reset();
	iterModel->reset();
	iterFlag->reset();
}
void VisAccum::nextTime()
{
	iterData->next();
	iterWeight->next();
	iterModel->next();
	iterFlag->next();	
}
bool VisAccum::accumulate(const Table& tab)
{
	Vector<Int> ant1,ant2,scan,spw,field;
	Vector<double> time;
	//double time;
	Cube<Complex> data,model;
	Cube<float> weight;
	Cube<bool> flag;
	Vector<int> bindx;
	Vector<Slice> tslice;
	ScalarColumn<Int> _ant1 (tab,  MS::columnName(MS::ANTENNA1));
	ScalarColumn<Int> _ant2 (tab, MS::columnName(MS::ANTENNA2));
	ScalarColumn<double> _time (tab, MS::columnName(MS::TIME));
	ScalarColumn<int> _scan (tab, MS::columnName(MS::SCAN_NUMBER));
	ScalarColumn<int> _spw (tab, MS::columnName(MS::DATA_DESC_ID));
	ScalarColumn<int> _field (tab, MS::columnName(MS::FIELD_ID));
	
	ArrayColumn<Complex> _data (tab,  MS::columnName(MS::DATA));
	
	ArrayColumn<Complex> _model;
	if(hasmodel)
		_model=ArrayColumn<Complex>(tab, MS::columnName(MS::MODEL_DATA));
	ArrayColumn<float> _weight;
	if(hasWtSpec)
	{
	 	_weight=ArrayColumn<float>(tab, MS::columnName(MS::SIGMA_SPECTRUM));
	}
	else
	{
		_weight=ArrayColumn<float>(tab, MS::columnName(MS::SIGMA));
	}
	ArrayColumn<bool> _flag (tab, MS::columnName(MS::FLAG));


	time.reference(_time.getColumn());
	scan.reference(_scan.getColumn());
	spw.reference(_spw.getColumn());
	field.reference(_field.getColumn());
	data.reference(_data.getColumn(cslice));
	//cout<<data.shape()<<endl;
	if(hasmodel)
		model.reference(_model.getColumn(cslice));
	flag.reference(_flag.getColumn(cslice));
	//cout<<hasWtSpec<<endl;
	if(hasWtSpec)
		weight.reference(_weight.getColumn(cslice));
	else
	{
		
		IPosition dim(3);
		dim[0]=data.shape()[0];
		dim[1]=1;
		dim[2]=data.shape()[2];
		Vector< Vector<Slice> > cslice_tmp;
		cslice_tmp.resize(1);
		cslice_tmp[0]=cslice[0];
		weight.reference(_weight.getColumn(cslice_tmp).reform(dim));
	}
	tslice=getTimeStrides(time);
	//cout<<tstart<<","<<tslice.shape()[0]<<endl;
	for (int i=tstart;i<tslice.shape()[0];i++)
	{	
		//cout<<"i:"<<i<<endl;
		//cout<<tslice[i]<<endl;
		bindx.reference(getBaselineIndx(_ant1.getColumnRange(tslice[i]),_ant2.getColumnRange(tslice[i])));

		accumTime+=time[tslice[i].start()];
		accumScan=scan[tslice[i].start()];
		accumSPW=spw[tslice[i].start()];
		accumField=field[tslice[i].start()];
		//cout<<accumScan<<","<<accumField<<endl;
		for(int k=0,irow=tslice[i].start();k<(bindx).shape()[0];k++,irow++)
		{
			//cout<<k<<","<<bindx[k]<<","<<irow<<endl;
			Cube<Complex>(iterData->array()).xyPlane((bindx)[k])=data.xyPlane(irow);	
			if(hasmodel)
				Cube<Complex>(iterModel->array()).xyPlane((bindx)[k])=model.xyPlane(irow);	
			Cube<float>(iterWeight->array()).xyPlane((bindx)[k])=weight.xyPlane(irow);
			Cube<bool>(iterFlag->array()).xyPlane((bindx)[k])=flag.xyPlane(irow);	
		}

		nsolint--;
		if(nsolint==0)
		{
			//cout<<accumSPW<<endl;
			tstart=i+1;
			return false;
		}
		nextTime();
	} 
	//cout<<Cube<float>(iterWeight->array()).shape()<<endl;
	//cout<<Cube<float>(iterWeight->array()).xyPlane(0)<<endl;
	tstart=0;
	return true;
}
void VisAccum::initIter(double timeInteval)
{
	Block<int> sort(1);
	//sort[1] = MS::SPECTRAL_WINDOW;
	//sort[0] = MS::ANTENNA1;
	//sort[1] = MS::ANTENNA2;
	sort[0] = MS::TIME;
	//sort[1] = MS::ARRAY_ID;
	//sort[2] = MS::DATA_DESC_ID;
	msIter=new MSIter(ms,sort,timeInteval,false,false);
	msIter->origin();
}
bool VisAccum::nextIter()
{
	int i=0;
	bool donext;	
	
	donext=accumulate(msIter->table());
	//cout<<"donext:"<<donext<<endl;
	
	if(donext)
		(*msIter)++;
	endflag=msIter->more();
	//cout<<"endflag:"<<endflag<<endl;
	return (!donext);
	
}
void VisAccum::setnsolint(int _nsolint)
{
	nsolint=_nsolint;
}
bool VisAccum::getEndflag()
{
	return endflag;
}
int VisAccum::getNant()
{
	return nant;
}
PyObject*  VisAccum::getAccum()
{
	PyObject *pyData,*pyModel,*pyWeight,*pyFlag;
	pyData=casacore::python::casa_array_to_python<Complex>::convert(accumData);
	pyModel=casacore::python::casa_array_to_python<Complex>::convert(accumModel);
	pyWeight=casacore::python::casa_array_to_python<float>::convert(accumWeight);
	pyFlag=casacore::python::casa_array_to_python<bool>::convert(accumFlag);
	return Py_BuildValue("diiiOOOO",accumTime,accumScan,accumSPW,accumField,pyData,pyModel,pyWeight,pyFlag);
}

BOOST_PYTHON_MODULE(VisAccum)
{
     class_<VisAccum>("VisAccum",init<string,string,string,string,string,string,string,int,double,int,int>())
        .def("resetAccum", &VisAccum::resetAccum)
        .def("nextIter", &VisAccum::nextIter)
	.def("setnsolint", &VisAccum::setnsolint)
	.def("getEndflag",&VisAccum::getEndflag)
	.def("getAccum",&VisAccum::getAccum)
	.def("getNant",&VisAccum::getNant)
	.def("setHasModel",&VisAccum::setHasModel)
    ;
}
