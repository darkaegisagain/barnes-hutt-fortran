barneshutt: barneshutt.f90
	ifort -qopenmp -c -m64 -O2 barneshutt.f90
	ifort -qopenmp -g -m64 -O2 -o barneshutt barneshutt.f90


