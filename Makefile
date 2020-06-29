FLAGS = -fPIC

# the python interface through swig
PYTHONI = -I/usr/include/python3.7m/
PYTHONL = -Xlinker -export-dynamic 
CFLAGS = -I/usr/include
LIBFLAGS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm
# default super-target
R: 
	g++ $(CFLAGS) $(LINFLAGS) -march=native -O3 -fPIC -c R.c -o R.o
	swig -c++ -python -o R_wrap.cxx R.i 
	g++ $(FLAGS) $(PYTHONI)  -c R_wrap.cxx -o R_wrap.o



muscle2pop: 
	g++ $(CFLAGS)  -march=native -O3 -fPIC -c $@.c -o $@.o
	swig -c++ -python -o $@_wrap.cxx $@.i 
	g++ $(FLAGS) $(PYTHONI)  -c $@_wrap.cxx -o $@_wrap.o
	g++ $@.o $@_wrap.o $(PYTHONL) $(LIBFLAGS) -shared  -o _$@.so

x2POP: 
	g++ $(CFLAGS)  -march=native -O3 -fPIC -c $@.c -o $@.o
	swig -c++ -python -o $@_wrap.cxx $@.i 
	g++ $(FLAGS) $(PYTHONI)  -c $@_wrap.cxx -o $@_wrap.o
	g++ $@.o $@_wrap.o $(PYTHONL) $(LIBFLAGS) -shared  -o _$@.so


x2POPF: 
	g++ $(CFLAGS)  -march=native -O3 -fPIC -c $@.c -o $@.o
	swig -c++ -python -o $@_wrap.cxx $@.i 
	g++ $(FLAGS) $(PYTHONI)  -c $@_wrap.cxx -o $@_wrap.o
	g++ $@.o $@_wrap.o $(PYTHONL) $(LIBFLAGS) -shared  -o _$@.so
