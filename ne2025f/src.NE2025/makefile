
OBJ = NE2025.o dmdsm.NE2025.o  density.NE2025.o neLISM.NE2025.o neclumpN.NE2025.o nevoidN.NE2025.o scattering98.o 
CFLAGS = 
BINDIR = ../bin.NE2025/


NE2025: $(OBJ)
	gfortran $(CFLAGS) $(OBJ) -o $(BINDIR)NE2025


NE2025.o: NE2025.f 

getarg.o: gerard.f

iargc.o: iargc.f

density.NE2025.o: density.NE2025.f

dmdsm.NE2025.o: dmdsm.NE2025.f

neLISM.NE2025.o: neLISM.NE2025.f

neclumpN.NE2025.o: neclumpN.NE2025.f

nevoidN.NE2025.o: nevoidN.NE2025.f


scattering98.o: scattering98.f

clean: 
	'rm' *.o
