#
#F90 = g95
F90 = gfortran -g 
#F90 = ifort

#newton: test_newton.o root_find.o protostellar_core.o constants.o
#	${F90} -o xtest_newton test_newton.o root_find.o protostellar_core.o constants.o
##
#test_newton.o: test_newton.f90 root_find.o protostellar_core.o constants.o
#	${F90} -c test_newton.f90
##
protostar: clean protostar.o rk2.o root_find.o protostellar_core.o constants.o
	${F90} -o xprotostar protostar.o rk2.o root_find.o protostellar_core.o constants.o
#
protostar.o: protostar.f90 rk2.o root_find.o protostellar_core.o constants.o
	${F90} -c protostar.f90
#
rk2.o: rk2.f90
	${F90} -c rk2.f90
#
root_find.o: root_find.f90 protostellar_core.o constants.o
	${F90} -c root_find.f90
#
protostellar_core.o: protostellar_core.f90 constants.o
	${F90} -c protostellar_core.f90
#
constants.o: constants.f90
	${F90} -c constants.f90


clean:
	-rm *.o *.mod xprotostar Beta_interpolation_file.txt
#
#machine_epsilon: test_machine_epsilon.o machine_epsilon.o constants.o
#	${F90} -o xtest_machine_epsilon test_machine_epsilon.o machine_epsilon.o constants.o
##
#test_machine_epsilon.o: test_machine_epsilon.f90 machine_epsilon.o constants.o
#	${F90} -c test_machine_epsilon.f90
##
#machine_epsilon.o: machine_epsilon.f90 constants.o
#	${F90} -c machine_epsilon.f90
##
##
