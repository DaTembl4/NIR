three: plot.plt two
	gnuplot -p plot.plt
two: one
	./a
one: period.f90 source.dat plot.plt
	gfortran period.f90


