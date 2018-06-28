TARGET = ene_lanczos_main_1114.f90 sz_z_0119.f90 sparse_mat_z_1114.f90 \
          lanczos_1114.f90 quicksort.f90 berry_phase_1114.f90 pair_format.f90 observable.f90
OBJECTS = ene_lanczos_main_1114.o sz_z_0119.o sparse_mat_z_1114.o \
          lanczos_1114.o quicksort.o berry_phase_1114.o pair_format.o  observable.o
#OBJECTS = ene_lanczos_z_0119.f90 sz_z_0119.f90 sparse_mat_z_0129_mp.f90 \
          lanczos_z_0119.f90 quicksort.f90 berry_phase_0123.f90
#OBJECTS = test.o
MOD_FILES  =
#FC = gfortran
FC = mpif90
#FC = ifort
# for gfortran
ifeq (${FC},gfortran)
	FFLAGS += -fimplicit-none -fbounds-check
	LDFLAGS += -fopenmp -llapack -lblas
endif
# for ifort
ifeq (${FC},ifort)
#	FFLAGS = -O2 -xCORE-AVX2 -qopenmp #-g -check all -traceback
	FFLAGS = -O2 # -g #-check all -traceback
	LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
endif
#for mpif90
ifeq (${FC},mpif90)
#	FFLAGS = -O2 -xCORE-AVX2 -qopenmp #-g -check all -traceback
	FFLAGS = -O2 -xCORE-AVX2 -qopenmp #-g #-check all -traceback
	LDFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
endif
test3:
	$(FC) $(FFLAGS) -c $(TARGET) # $(LDFLAGS)
	$(FC) $(FFLAGS) -o test3 $(OBJECTS) $(LDFLAGS)
spin2-1:
	$(FC) $(FFLAGS) -c $(TARGET) # $(LDFLAGS)
	$(FC) $(FFLAGS) -o spin2-1 $(OBJECTS) $(LDFLAGS)
spin2-2:
	$(FC) $(FFLAGS) -c $(TARGET) # $(LDFLAGS)
	$(FC) $(FFLAGS) -o spin2-2 $(OBJECTS) $(LDFLAGS)
spin2-3:
	$(FC) $(FFLAGS) -c $(TARGET) # $(LDFLAGS)
	$(FC) $(FFLAGS) -o spin2-3 $(OBJECTS) $(LDFLAGS)
spin2-4:
	$(FC) $(FFLAGS) -c $(TARGET) # $(LDFLAGS)
	$(FC) $(FFLAGS) -o spin2-4 $(OBJECTS) $(LDFLAGS)
spin2-5:
	$(FC) $(FFLAGS) -c $(TARGET) # $(LDFLAGS)
	$(FC) $(FFLAGS) -o spin2-5 $(OBJECTS) $(LDFLAGS)
spin2-6:
	$(FC) $(FFLAGS) -c $(TARGET) # $(LDFLAGS)
	$(FC) $(FFLAGS) -o spin2-6 $(OBJECTS) $(LDFLAGS)

# .PHONY: clean
 clean:
	${RM} ${OBJECTS} ${MOD_FILES} test3 spin2-*
