NAME=PC
#NAME=LAPTOP
#NAME=RCAC
ifeq ($(NAME), RCAC)   # Hansen
INC=/home/cheng109/work/phosim_DES/phosim/phosim_core/source/cfitsio/include
LIB=/home/cheng109/work/phosim_DES/phosim/phosim_core/source/cfitsio/lib
INC_EIGEN=/home/cheng109/share/eigen_include/
endif
ifeq ($(NAME), PC)
INC=-I/usr/local/Cellar/boost/1.60.0_1/include -I/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/include -I/usr/local/include 
LIB=-L/usr/local/Cellar/boost/1.60.0_1/lib -L/usr/local/lib -L/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/lib -L/opt/local/lib -L/usr/local/gfortran/lib
INC_EIGEN=-I/Users/cheng109/eigen_include
endif
ifeq ($(NAME), LAPTOP)
INC=-I/usr/local/Cellar/boost/1.58.0/include -I/Users/juncheng/work/cfitsio -I/usr/local/include -I/usr/local/Cellar/gsl/1.16/include
LIB=-L/usr/local/Cellar/boost/1.58.0/lib -L/Users/juncheng/work/cfitsio -L/opt/local/lib -L/usr/local/Cellar/gsl/1.16/lib -L/usr/local/gfortran/lib
#-L/usr/local/Cellar/gcc/5.2.0/lib/gcc/5
INC_EIGEN=-I/Users/juncheng/work/eigen_include
endif
CC=g++
#gcc
#clang++
CFLAGS=-Wall -O3 $(INC) $(INC_EIGEN) $(LIB)
LDFLAGS=-lcfitsio -lgsl -L. -lfortranstuff -lgfortran 
#-lg2c
# -larmadillo -lboost_iostreams -lboost_system #-fopenmp 
SRCS=main.cpp Image.cpp commons.cpp Model.cpp
FC=gfortran

COMMON_HDRS= Image.h Model.h commons.h fastell.h parafit.h  
OBJS= commons.o main.o Image.o Model.o parafit.o
#OMP_NUM_THREADS=8
#F77FLAGS := -g $(OPT) -c  -ffixed-line-length-none -fno-automatic -fbounds-check 
#-fno-underscoring
#all: main.o Image.o commons.o Model.o
#	g++ $(INC) $(LIB) $(CFLAG) $(LDFLAGS) commons.cpp main.cpp Image.cpp Model.cpp -o junGL 
#	./junGL


## ARG1:   working directory. 
## ARG2:   configuration file. 
## ARG3:   output files. 
ARG1=horseshoe_test/
ARG2=conf.txt
ARG3=output.txt

all:  $(COMMON_HDRS) $(OBJS) libfortranstuff.a

	@#$(CC) $(CFLAGS) -o junGL commons.cpp main.cpp Image.cpp Model.cpp gl_crit.cpp parafit.cpp $(LDFLAGS)
	@$(CC) $(CFLAGS) -o junGL $(OBJS) $(LDFLAGS)
	@#valgrind --tool=memcheck --leak-check=full --verbose --log-file=memcheck.log --track-origins=yes ./junGL
	@./junGL $(ARG1) $(ARG2) $(ARG3)
commons.o: commons.cpp
	@$(CC) -c -o commons.o commons.cpp $(CFLAGS)
main.o:  main.cpp
	@$(CC) -c -o main.o main.cpp $(CFLAGS)
Image.o:  Image.cpp Image.h
	@$(CC) -c -o Image.o Image.cpp $(CFLAGS)
Model.o: Model.cpp
	@$(CC) -c -o Model.o Model.cpp $(CFLAGS)
parafit.o: parafit.cpp
	@$(CC) -c -o parafit.o parafit.cpp $(CFLAGS)
libfortranstuff.a:
	@$(FC) -O -c slatec/src/*.f fastell.f 
	@ar -r libfortranstuff.a *.o
	@#rm *.o
plot:
	@./plotScript
clean: 
	rm -f *.o junGL output.txt
imfit: 
	cd imfit-1.3; 
	scons imfit-1.3/imfit; 
#test:
#	$(FC) slatec/src/*.f fastell_example.f -o testFastEll
#	./testFastEll