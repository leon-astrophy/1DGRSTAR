include ./make.inc

SOURCES=Definition.f90 \
	Main.f90 \
	BuildHydro.f90 \
	EOSTABLE.f90 \
	EOS-Analytic.f90 \
	FermiMo.f90 \
	GetRho.f90 \
	GetRho_EOSPToR.f90 \
	GetRho_EOSRToP.f90 \
	GetEpsilon.f90 \
	Ini_Der.f90 \
	NSEOS.f90 

MODULES=Definition.mod
OBJECTS=$(SOURCES:.f90=.o )

TOV: Definition.o $(OBJECTS)  
	$(F90) $(LDFLAGS) -o ./TOV $(OBJECTS) 

$(OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) -c $< -o $@

TOV.o: Definition.o

clean:
	rm -rf Definition
	rm -rf *.o 	
	rm -rf *.mod
cleanfile:
	rm -rf ./Outfile/*.dat
	rm -rf TOV
	rm -rf tmp.txt
