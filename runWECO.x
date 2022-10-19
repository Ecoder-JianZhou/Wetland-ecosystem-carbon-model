#/bin/sh
gfortran -g src/dataTypes.F90 src/driver.f90 src/vegetation.f90 src/soil.f90 src/transferCpools.f90 src/methane.f90 src/main.f90 -o weco

./weco

rm weco

