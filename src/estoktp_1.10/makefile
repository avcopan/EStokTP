#FLIBS = -lblas
estoktp.x: estoktp.o g09fopt.o LineRead.o  mt.o croutines.o molprofopt.o zmat_routines.o
	gfortran  -o ../../exe/estoktp.x  estoktp.o g09fopt.o LineRead.o croutines.o mt.o molprofopt.o zmat_routines.o
estoktp.o: estoktp.f
	gfortran -c  estoktp.f -o estoktp.o 
g09fopt.o: g09fopt.f
	gfortran  -c g09fopt.f -o g09fopt.o 
LineRead.o: LineRead.f
	gfortran  -c LineRead.f -o LineRead.o 
molprofopt.o: molprofopt.f
	gfortran -c  molprofopt.f -o molprofopt.o 
mt.o: mt.c
	gcc -c -fno-leading-underscore mt.c -o mt.o  
croutines.o: croutines.c
	gcc -c -fno-leading-underscore  croutines.c -o croutines.o  
zmat_routines.o: zmat_routines.f
	gfortran -c zmat_routines.f -o zmat_routines.o

