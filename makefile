all:SimGoldbeter
SimGoldbeter:	
	gfortran -o 3CompCellMov main-Source.f anfang.f rs-KeEverywhere.f FourPillars.f \
	out.f ODE-Merson.f ic-Source.f StartingTime.f LapAndMovement.f \
	-L/usr/bmp/pgplot-5.2/ -lpgplot \
	-L/usr/bmp/slatec-4.1/lib -lslatec \
	-L/usr/bmp/lapack-3.4.0 \
	-lX11 \
	-lpng
