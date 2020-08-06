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
using namespace casacore;
using namespace boost::python;
class VisAccum
{
	public:
		VisAccum(string msname,string fieldname, string spw,string uvrange,string scan,string observation,string poln,int ntimesample,double _timeInverval,int dataDescIndx,int spwIndx);
		void resetAccum();
		void setnsolint(int nsolint);
		bool getEndflag();
		bool nextIter();
		int getNant();
		void setHasModel(bool _hasmodel);
		PyObject*  getAccum();
		
	private:
		int nsolint;
		IPosition getDim(int ntimesample);
		IPosition dim;
		void msSelect(string &msname,string &fieldname, string &spw,string &uvrange,string &scan,string &observation,string& poln,int dataDescIndx,int spwIndx);
		Vector<int> getBaselineIndx(const Vector<int>& ant1,const Vector<int>& ant2);
		Vector<Slice> getTimeStrides(const Vector<double>& time);
		void initAccum(IPosition dim);
		void nextTime();
		void initIter(double timeInteval);
		bool accumulate(const Table& tab);
		bool hasmodel;
		bool hasWtSpec;
		MeasurementSet ms;
		Array<Complex> accumData,accumModel;
		Array<float> accumWeight;
		Array<bool> accumFlag;
		MSIter* msIter;
		double accumTime;
		bool endflag;
		int tstart;
		int accumScan;
		int accumSPW,accumField;
		Vector<int> solintMap;
		ArrayIterator<Complex>* iterData;
		ArrayIterator<float>* iterWeight;
		ArrayIterator<Complex>* iterModel;
		ArrayIterator<bool>* iterFlag;
		Vector< Vector<Slice> > cslice;
		int nant;
	
};
