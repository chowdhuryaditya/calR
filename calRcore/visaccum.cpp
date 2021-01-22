#include<visaccum.h>
using namespace casacore;
using namespace boost::python;
VisAccum::VisAccum(string msname,string fieldname, string spw,string uvrange,string scan,string observation,string poln,int ntimesample,double timeInverval,int dataDescIndx,int spwIndx)
{
	msSelect(msname,fieldname,spw,uvrange,scan,observation,poln,dataDescIndx,spwIndx);
	MeasurementSet msOriginal(msname);
	nant=ROMSAntennaColumns(msOriginal.antenna()).nrow();
	nrowstepMax=10000;
	tstart=0;
	curRow=0;
	nrow=ms.nrow();
	//cout<<"nrow:"<<nrow<<endl;
	if(nrow<nrowstepMax)
		nrowstepMax=nrow;
	initAccum(getDim(ntimesample));	
	resetAccum();

	endflag=1;
	endreadflag=1;
	hasmodel=false;
	
}
void VisAccum::setHasModel(bool _hasmodel)
{
	hasmodel=_hasmodel;
	initIter();	
}
IPosition VisAccum::getDim(int ntimesample)
{
	Block<int> sort(1);
	sort[0] = MS::TIME;
	IPosition dim(4);
	dim[2]=(nant*(nant-1))/2;
	
	dim[3]=ntimesample;
	ArrayColumn<Complex> _data (ms,  MS::columnName(MS::DATA));
	Array<Complex> data=_data.getSlice(0,cslice);
	//cout<<data.shape()<<endl;
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
			tslice[itime]=Slice(curRow+startindxt,i-startindxt,1);
			curt=time[i];
			itime++;
			startindxt=i;
		}
		i++;
		
	}
	tslice[itime]=Slice(curRow+startindxt,i-startindxt,1);
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
	 	_weight=ArrayColumn<float>(ms, MS::columnName(MS::SIGMA_SPECTRUM));
		IPosition tmp=_weight.shape(0); //CHECK IMPLICATIONS OF THIS LINE
	}
	catch(...)
	{
		hasWtSpec=false;
		_weight=ArrayColumn<float>(ms, MS::columnName(MS::SIGMA));
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
long int VisAccum::getnRow()
{

	Vector<double> time2;
	ScalarColumn<double> _time (ms, MS::columnName(MS::TIME));
	Slicer rsel=Slicer(Slice(curRow,nrowstepMax));
	time2.reference(_time.getColumnRange(rsel));
	long int length=time2.shape()[0];
	for(long int i=length-1;i>0;i--)
	{
		if(time2[i]!=time2[i-1])
		{
			return i;
		}
	}
	return 0;
	
}
void VisAccum::readFromMS()
{
	//cout<<nrowstep<<endl;
	data.reference(Cube<Complex>(accumData.shape()[0],accumData.shape()[1],nrowstep));
	model.reference(Cube<Complex>(accumData.shape()[0],accumData.shape()[1],nrowstep));
	flag.reference(Cube<bool> (accumData.shape()[0],accumData.shape()[1],nrowstep));
	if(hasWtSpec)
		weight.reference(Cube<float>(accumData.shape()[0],accumData.shape()[1],nrowstep));
	else
		weight.reference(Cube<float>(accumData.shape()[0],1,nrowstep));

	//cout<<data.shape()<<endl;
	ScalarColumn<Int> _ant1 (ms,  MS::columnName(MS::ANTENNA1));
	ScalarColumn<Int> _ant2 (ms, MS::columnName(MS::ANTENNA2));
	ScalarColumn<double> _time (ms, MS::columnName(MS::TIME));
	ScalarColumn<int> _scan (ms, MS::columnName(MS::SCAN_NUMBER));
	ScalarColumn<int> _spw (ms, MS::columnName(MS::DATA_DESC_ID));
	ScalarColumn<int> _field (ms, MS::columnName(MS::FIELD_ID));
	
	ArrayColumn<Complex> _data (ms,  MS::columnName(MS::DATA));
	
	ArrayColumn<Complex> _model;
	
	if(hasmodel)
		_model=ArrayColumn<Complex>(ms, MS::columnName(MS::MODEL_DATA));
	ArrayColumn<float> _weight;
	if(hasWtSpec)
	{
	 	_weight=ArrayColumn<float>(ms, MS::columnName(MS::SIGMA_SPECTRUM));
	}
	else
	{
		_weight=ArrayColumn<float>(ms, MS::columnName(MS::SIGMA));
	}
	ArrayColumn<bool> _flag (ms, MS::columnName(MS::FLAG));


	time.reference(_time.getColumnRange(rowselect));

	scan.reference(_scan.getColumnRange(rowselect));
	spw.reference(_spw.getColumnRange(rowselect));
	field.reference(_field.getColumnRange(rowselect));
	
	for(int irow=0;irow<nrowstep;irow++)
	{
		data.xyPlane(irow)=_data.getSlice(irow+curRow,cslice);
		flag.xyPlane(irow)=_flag.getSlice(irow+curRow,cslice);
		
	}
	if(hasmodel)
	{
		for(int irow=0;irow<nrowstep;irow++)
		{
			model.xyPlane(irow)=_model.getSlice(irow+curRow,cslice);
		}
	}
	if(hasWtSpec)
	{
		for(int irow=0;irow<nrowstep;irow++)
		{
			weight.xyPlane(irow)=_weight.getSlice(irow+curRow,cslice);
		}
	}
	else
	{
			
		IPosition dim(2);
		dim[0]=data.shape()[0];
		dim[1]=1;
		Vector< Vector<Slice> > cslice_tmp;
		cslice_tmp.resize(1);
		cslice_tmp[0]=cslice[0];
		for(int irow=0;irow<nrowstep;irow++)
		{
			weight.xyPlane(irow)=_weight.getSlice(irow+curRow,cslice_tmp).reform(dim);
		}
	}
}


bool VisAccum::accumulate()
{

	
	//cout<<hasWtSpec<<endl;

	Vector<int> bindx;
	Vector<Slice> tslice;
	ScalarColumn<Int> _ant1 (ms,  MS::columnName(MS::ANTENNA1));
	ScalarColumn<Int> _ant2 (ms, MS::columnName(MS::ANTENNA2));
	tslice=getTimeStrides(time);
	int sindx;
	//cout<<tstart<<","<<tslice.shape()[0]<<endl;
	for (int i=tstart;i<tslice.shape()[0];i++)
	{	
		//cout<<"i:"<<i<<endl;
		//cout<<tslice[i]<<endl;
		bindx.reference(getBaselineIndx(_ant1.getColumnRange(tslice[i]),_ant2.getColumnRange(tslice[i])));
		sindx=tslice[i].start()-curRow;
		accumTime+=time[sindx];
		accumScan=scan[sindx];
		accumSPW=spw[sindx];
		accumField=field[sindx];
		//cout<<accumScan<<","<<accumField<<endl;
		for(int k=0,irow=sindx;k<(bindx).shape()[0];k++,irow++)
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
void VisAccum::initIter()
{
	//Block<int> sort(1);
	//sort[1] = MS::SPECTRAL_WINDOW;
	//sort[0] = MS::ANTENNA1;
	//sort[1] = MS::ANTENNA2;
	//sort[0] = MS::TIME;
	//sort[1] = MS::ARRAY_ID;
	//sort[2] = MS::DATA_DESC_ID;
	//msIter=new MSIter(ms,sort,timeInteval,false,false);
	//msIter->origin();
	curRow=0;
	nrowstep=getnRow();
	rowselect=Slicer(Slice(0,nrowstep));
	readFromMS();
	
}
bool VisAccum::nextIter()
{
	int i=0;
	bool donext;	
	
	donext=accumulate();
	//cout<<"donext:"<<donext<<endl;
	if(donext and !endreadflag)
	{
		endflag=0;
	}
	if(donext and endreadflag)
	{
		curRow+=nrowstep;
		if(curRow+nrowstepMax>nrow)
		{
			nrowstep=nrow-curRow;
			endreadflag=0;
		}
		else
		{
			nrowstep=getnRow();
		}
		rowselect=Slicer(Slice(curRow,nrowstep));
		readFromMS();

	}
	
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
boost::python::object VisAccum::getData()
{
	PyObject *pyData,*pyModel,*pyWeight,*pyFlag;
	pyData=casacore::python::casa_array_to_python<Complex>::convert(accumData);
	
	boost::python::handle<> handleData(pyData);
	
	return boost::python::object(handleData);
}
boost::python::object VisAccum::getModel()
{
	PyObject *pyModel;
	
	pyModel=casacore::python::casa_array_to_python<Complex>::convert(accumModel);
	boost::python::handle<> handleModel(pyModel);
	
	return boost::python::object(handleModel);
}
boost::python::object VisAccum::getWeight()
{
	PyObject *pyWeight;
	pyWeight=casacore::python::casa_array_to_python<float>::convert(accumWeight);
	boost::python::handle<> handleWeight(pyWeight);
	return boost::python::object(handleWeight);
}
boost::python::object VisAccum::getFlag()
{
	PyObject *pyFlag;
	pyFlag=casacore::python::casa_array_to_python<bool>::convert(accumFlag);	
	boost::python::handle<> handleFlag(pyFlag);
	return boost::python::object(handleFlag);
}
PyObject*  VisAccum::getAccum()
{
	PyObject* PyReturnval= Py_BuildValue("diii",accumTime,accumScan,accumSPW,accumField);
	return PyReturnval;
	
}
BOOST_PYTHON_MODULE(VisAccum)
{
     class_<VisAccum>("VisAccum",init<string,string,string,string,string,string,string,int,double,int,int>())
        .def("resetAccum", &VisAccum::resetAccum)
        .def("nextIter", &VisAccum::nextIter)
	.def("setnsolint", &VisAccum::setnsolint)
	.def("getEndflag",&VisAccum::getEndflag)

	.def("getAccum",&VisAccum::getAccum)
	.def("getData",&VisAccum::getData)
	.def("getModel",&VisAccum::getModel)
	.def("getWeight",&VisAccum::getWeight)
	.def("getFlag",&VisAccum::getFlag)

	.def("getNant",&VisAccum::getNant)
	.def("setHasModel",&VisAccum::setHasModel)
    ;
}

