# Binary file
BIN			= 2eclvsB

# Flags
CFLAGS			= -O3 -cpp
LDFLAGS			= -llapack #-mkl=sequential -lmkl_lapack95_lp64 -lpthread -lm

# Default Values
zmin_			= -2000.d0
zmax_			= 2000.d0
z0_			= 0.d0
sigma_			= 20.d0
lum_			= 0
c_			= 0.d0
gamma_			= 0.d0
kord_			= 5
l_			= 70
intg_			= 60
nev_			= 50
V0_			= 0.05d0
B_campo_i_		= 10.d0
B_campo_f_		= 100.d0
num_puntos_B_		= 90
me_			= 0.063d0
lambda_			= 0.4d0

# OFILE			= [CPU,$(Q),$(L),$(SAMPLES),$(TEMP_MIN),$(TEMP_MAX),$(DELTA_TEMP),$(TRAN),$(TMAX),$(DELTA_T)].dat

# Simulation Parameters
PARAMETERS		= -Dzmin_=$(zmin_) -Dzmax_=$(zmax_) -Dz0_=$(z0_) -Dsigma_=$(sigma_) -Dlum_=$(lum_) -Dc_=$(c_) \
			  -Dgamma_=$(gamma_) -Dkord_=$(kord_) -Dl_=$(l_) -Dintg_=$(intg_) -Dnev_=$(nev_) -DV0_=$(V0_) \
			  -DB_campo_i_=$(B_campo_i_) -DB_campo_f_=$(B_campo_f_) -Dnum_puntos_B_=$(num_puntos_B_) \
			  -Dme_=$(me_) -Dlambda_=$(lambda_)

# Compilers
CC			= gfortran #ifort
LINKER			= gfortran #ifort

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
