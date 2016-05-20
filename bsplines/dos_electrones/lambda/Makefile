# Binary file
BIN			= 2eclvslambda

# Flags
# CFLAGS			= -O3 -Wall -Wextra -cpp #-fopenmp
CFLAGS			= -O3 -cpp
LDFLAGS			= -llapack #-mkl=sequential -lmkl_lapack95_lp64 -lpthread -lm

# Default Values
zmin_			= -100.d0
zmax_			= 100.d0
z0_			= 0.d0
sigma_			= 20.d0
lum_			= 0
c_			= 0.d0
gamma_			= 0.d0
kord_			= 5
l_			= 60
intg_			= 50
nev_			= 50
V0_			= 0.05d0
lambdai_		= 0.d0
lambdaf_		= 1.d0
num_puntos_lambda_	= 20
me_			= 0.063d0
B_campo_		= 80.d0

# OFILE			= [CPU,$(Q),$(L),$(SAMPLES),$(TEMP_MIN),$(TEMP_MAX),$(DELTA_TEMP),$(TRAN),$(TMAX),$(DELTA_T)].dat

# Simulation Parameters
PARAMETERS		= -Dzmin_=$(zmin_) -Dzmax_=$(zmax_) -Dz0_=$(z0_) -Dsigma_=$(sigma_) -Dlum_=$(lum_) -Dc_=$(c_) \
			  -Dgamma_=$(gamma_) -Dkord_=$(kord_) -Dl_=$(l_) -Dintg_=$(intg_) -Dnev_=$(nev_) -DV0_=$(V0_) \
			  -Dlambdai_=$(lambdai_) -Dlambdaf_=$(lambdaf_) -Dnum_puntos_lambda_=$(num_puntos_lambda_) -Dme_=$(me_) -DB_campo_=$(B_campo_)

# Compilers
CC			= gfortran-4.8 #ifort
LINKER			= gfortran-4.8 #ifort

# Files
C_SOURCES		= $(BIN).f90
HEADERS			=
C_OBJS			= $(patsubst %.f90, %.o, $(C_SOURCES))

# Rules

$(BIN): clean $(C_OBJS) $(HEADERS)
	$(LINKER) $(CFLAGS) -o $(BIN) $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CC) -c $(C_SOURCES) $(CFLAGS) $(INCLUDES) $(PARAMETERS)

run: $(MP) $(BIN)
	./$(BIN) #> $(OFILE) &

clean:
	rm -f $(BIN) *.o *.mod
