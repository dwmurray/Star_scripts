#
#F90 = g95
F90 = gfortran
#F90 = ifort
machine_epsilon: test_machine_epsilon.o machine_epsilon.o constants.o
	${F90} -o xtest_machine_epsilon test_machine_epsilon.o machine_epsilon.o constants.o
#
test_machine_epsilon.o: test_machine_epsilon.f90 machine_epsilon.o constants.o
	${F90} -c test_machine_epsilon.f90
#
machine_epsilon.o: machine_epsilon.f90 constants.o
	${F90} -c machine_epsilon.f90
#
#
constants.o: constants.f90
	${F90} -c constants.f90
