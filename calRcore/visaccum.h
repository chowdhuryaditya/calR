#include<string>
#include <boost/python.hpp>
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSSelectionTools.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSIter.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/python/Converters/PycArray.h>
using namespace std;

class VisAccum
{
	public:
		VisAccum(string msname,string fieldname, string spw,string uvrange,string scan,string observation,string poln,int ntimesample,double timeInverval,int dataDescIndx,int spwIndx);
		void resetAccum();
		void setnsolint(int nsolint);
		bool getEndflag();
		bool nextIter();
		int getNant();
		void setHasModel(bool _hasmodel);
		boost::python::object getData();
		boost::python::object getModel();
		boost::python::object getWeight();
		boost::python::object getFlag();
		PyObject*  getAccum();
		
	private:
		casacore::Vector<casacore::Int> ant1,ant2,scan,spw,field;
		casacore::Vector<double> time;
		casacore::Cube<casacore::Complex> data;
		casacore::Cube<casacore::Complex> model;
		casacore::Cube<bool> flag;
		casacore::Cube<float> weight;

		int nsolint;
		unsigned long int curRow;
		unsigned long int nrow;
		unsigned int nrowstepMax;
		unsigned int nrowstep;
		casacore::Slicer rowselect;
		double timeInterval;
		casacore::IPosition getDim(int ntimesample);
		casacore::IPosition dim;
		void msSelect(string &msname,string &fieldname, string &spw,string &uvrange,string &scan,string &observation,string &poln,int dataDescIndx,int spwIndx);
		casacore::Vector<int> getBaselineIndx(const casacore::Vector<int>& ant1,const casacore::Vector<int>& ant2);
		casacore::Vector<casacore::Slice> getTimeStrides(const casacore::Vector<double>& time);
		void initAccum(casacore::IPosition dim);
		void nextTime();
		void readFromMS();
		void initIter();
		long int getnRow();
		bool accumulate();
		bool hasmodel;
		bool hasWtSpec;
		casacore::MeasurementSet ms;
		casacore::Array<casacore::Complex> accumData,accumModel;
		casacore::Array<float> accumWeight;
		casacore::Array<bool> accumFlag;
		//casacore::MSIter* msIter;
		double accumTime;
		bool endflag,endreadflag;
		int tstart;
		int accumScan;
		int accumSPW,accumField;
		casacore::Vector<int> solintMap;
		casacore::ArrayIterator<casacore::Complex>* iterData;
		casacore::ArrayIterator<float>* iterWeight;
		casacore::ArrayIterator<casacore::Complex>* iterModel;
		casacore::ArrayIterator<bool>* iterFlag;
		casacore::Vector< casacore::Vector<casacore::Slice> > cslice;
		casacore::Slicer cslicer;
		int nant;
	
};
