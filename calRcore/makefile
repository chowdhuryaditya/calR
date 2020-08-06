all:
	g++ -O3 -shared  -fPIC solver.cpp -o solver.so -L/usr/lib64 -lm -lstdc++ #-static  
	gcc -O3 -shared    -fPIC  visaccum.cpp -o VisAccum.so  -I. -I/usr/include/python2.7 -I/Data/chowdhury/soft/include libboost_python-mt.so libcasa_casa.so.2 libcasa_tables.so.2 libcasa_measures.so.2 libcasa_python.so.2 libcasa_ms.so.2 libcasa_meas.so.2 libcasa_scimath.so.2   -lstdc++ -lm  
