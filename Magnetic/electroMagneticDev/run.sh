gfortran -c externalControl.f -o control.o
gfortran -c meshElectroMagnetic.f -o mesh.o
gfortran -c fieldPhysical.f -o field.o
gfortran -c materialElectroMagnetic.f -o material.o
gfortran -c solverCalculate.f -o solver.o
gfortran -c postElectroMagnetic.f -o postElec.o
gfortran -g control.o mesh.o field.o material.o solver.o postElec.o main.f -o main
