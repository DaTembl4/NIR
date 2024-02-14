four: one two three
	./SKATSe
three: one two
	gfortran period.o uneventoeven.o -o SKATSe
two: uneventoeven.f90
	gfortran -c uneventoeven.f90 
one: period.f90
	gfortran -c period.f90


