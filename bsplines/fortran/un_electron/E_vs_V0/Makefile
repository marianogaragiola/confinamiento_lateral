# Binary file
BIN			= electron_Vlong

# Flags
CFLAGS			= -O3 -Wall -Wextra -cpp
LDFLAGS			= -llapack

# Default Values
zmin_			= -80.d0
zmax_			= 80.d0
z0_			= 0.d0
sigma_			= 15.d0
lum_			= 0
c_			= 0.d0
gamma_			= 0.d0
kord_			= 5
l_			= 90
intg_			= 100
nev_			= 15
me_			= 0.063d0
V0i_			= 0.d0
V0f_			= 0.1d0
num_puntos_V0_		= 100

# OFILE			= [CPU,$(Q),$(L),$(SAMPLES),$(TEMP_MIN),$(TEMP_MAX),$(DELTA_TEMP),$(TRAN),$(TMAX),$(DELTA_T)].dat

# Simulation Parameters
PARAMETERS		= -Dzmin_=$(zmin_) -Dzmax_=$(zmax_) -Dz0_=$(z0_) -Dsigma_=$(sigma_) -Dlum_=$(lum_) -Dc_=$(c_) \
			  -Dgamma_=$(gamma_) -Dkord_=$(kord_) -Dl_=$(l_) -Dintg_=$(intg_) -Dnev_=$(nev_) -Dme_=$(me_) \
			  -DV0i_=$(V0i_) -DV0f_=$(V0f_) -Dnum_puntos_V0_=$(num_puntos_V0_)

# Compilers
CC			= gfortran
LINKER			= gfortran

# Files
C_SOURCES		= $(BIN).f90
HEADERS			=
C_OBJS			= $(patsubst %.f90, %.o, $(C_SOURCES))


# Rules
$(BIN): clean $(C_OBJS) $(HEADERS)
	$(LINKER) $(CFLAGS) -o $(BIN) $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CC) -c $(C_SOURCES) $(CFLAGS) $(INCLUDES) $(PARAMETERS)

run: $(BIN)
	./$(BIN) #> $(OFILE) &

clean:
	rm -f $(BIN) *.o *.mod
