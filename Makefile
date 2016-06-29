OBJS = \
	subs.o\
	#sfuns.o
FFILE=set3d
#FLAGS=-O3 -fdefault-real-8
#FLAGS=-Wall -fcheck='all' -ffixed-line-length-none -fdefault-real-8 -O3 -J./Lib_VTK_IO/mod/
#FLAGS=-Wall -Wextra -Wconversion -fcheck='all' -ffixed-line-length-none -fdefault-real-8 -O3 -J ./Lib_VTK_IO/mod/
FLAGS=-ffixed-line-length-none -fdefault-real-8 -J ./Lib_VTK_IO/mod/
#STLFILE=twoCube10.stl
STLFILE=cube40.stl
#STLFILE=Ellipsoid2_coarse.stl
VTUFILE=test.vtu

all: compile run paraview

compile: VTK $(OBJS)
	#gfortran ${FFILE}.f90 ${OBJS} $(FLAGS) -o $(FFILE).exec
	gfortran $(FFILE).f90 ./Lib_VTK_IO/obj/*.o *.o $(FLAGS) -o $(FFILE).exec

file: 
	vim -O Makefile ../LevelSet_SA/Makefile

run:
	./$(FFILE).exec $(STLFILE)

VTK:
	cd Lib_VTK_IO; make

paraview:
	paraview $(VTUFILE) &

.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	gfortran $< $(FLAGS) -c -o $@

clean: 
	rm -f *.mod
	rm -f *.exec 
	rm -f *.o
	rm -f *.txt
	rm -rf output
