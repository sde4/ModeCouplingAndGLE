# Makefile for parallel tests     
# Lines starting with "#" are comments.
# To compile this parallel test, just type "make".

CC		= gcc
OMP		= -fopenmp
CFLAGS		= -O3 -g
LIBS		= -lgsl -lgslcblas -lm -lfftw3

# Stampede
FFT_PATH	= 
FFT_LIB  	= -L/opt/apps/gcc7_1/impi17_0/fftw3/3.3.6/lib
GSL_LIB  	= -L/opt/apps/gcc7_1/gsl/2.3/lib
FFT_INC		= -I/opt/apps/gcc7_1/impi17_0/fftw3/3.3.6/include
GSL_INC  	= -I/opt/apps/gcc7_1/gsl/2.3/include


# Golub
#FFT_PATH	= 
#FFT_LIB  	= -L/home/sde4/Programs/packages/fftw-3.3.1/lib
#GSL_LIB  	= -L/usr/local/gsl/2.4/lib
#FFT_INC		= -I/home/sde4/Programs/packages/fftw-3.3.1/include
#GSL_INC  	= -I/usr/local/gsl/2.4/include

SRC		= main_ModeCouplingAndGLE.c SysInit.c SysRead.c StateInit.c StateRead.c ModeCombinations4NonZeroCouplingConstants.c ModeCombinations4InternalResonance.c ThirdOrderCouplingCalculation.c ModeshapesDisp.c GalerkinProj.c IntegrateSys.c ForceSys.c StateStore.c array_def.c

OBJ		= $(SRC:.c=.o)

#EXEC1		= mcgle.serial
EXEC2		= mcgle.omp

all: $(SRC)  $(EXEC1)  $(EXEC2)

#$(EXEC1): $(OBJ)
#	$(CC) $(FFT_LIB) $(GSL_LIB) $(OBJ) $(LIBS) -o $(EXEC1)

$(EXEC2): $(OBJ)
	$(CC) $(OMP) $(FFT_LIB) $(GSL_LIB) $(OBJ) $(LIBS)  -o $(EXEC2)

%.o:%.c
	$(CC) $(OMP) $(CFLAGS)  $(GSL_INC) $(FFT_INC) -c $<
clean:
	rm -rf $(EXEC1) $(EXEC2)  *~ *.o core

cleanmode:
	rm -rf *.txt
cleaneng:
	rm -rf modeteng*
