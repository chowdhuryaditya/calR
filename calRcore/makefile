all:
	gcc -O3 -shared -fPIC solver.cpp -o solver.so -lm -lstdc++ 
	gcc -O3 -shared    -fPIC  visaccum.cpp -o VisAccum.so  -I. libboost_python-mt.so libcasa_casa.so.2 libcasa_tables.so.2 libcasa_measures.so.2 libcasa_python.so.2 libcasa_ms.so.2 libcasa_meas.so.2 libcasa_scimath.so.2 -lm -lstdc++ 

