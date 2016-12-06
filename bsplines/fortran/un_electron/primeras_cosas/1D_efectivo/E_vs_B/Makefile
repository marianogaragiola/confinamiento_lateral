# Binary file
BIN			= electron_B

# Flags
CFLAGS			= -O3 -Wall -Wextra -cpp
LDFLAGS			= -llapack

# Default Values
zmin_			= -100.d0
zmax_			= 100.d0
z0_			= 0.d0
sigma_			= 20.d0
lum_			= 0
c_			= 0.d0
gamma_			= 0.d0
kord_			= 5
l_			= 100
intg_			= 1000
nev_			= 30
me_			= 0.063d0
V0_			= 0.05d0
B_i_			= 5.d0
B_f_			= 100.d0
num_puntos_B_		= 190

# OFILE			= [CPU,$(Q),$(L),$(SAMPLES),$(TEMP_MIN),$(TEMP_MAX),$(DELTA_TEMP),$(TRAN),$(TMAX),$(DELTA_T)].dat

# Simulation Parameters
PARAMETERS		= -Dzmin_=$(zmin_) -Dzmax_=$(zmax_) -Dz0_=$(z0_) -Dsigma_=$(sigma_) -Dlum_=$(lum_) -Dc_=$(c_) \
			  -Dgamma_=$(gamma_) -Dkord_=$(kord_) -Dl_=$(l_) -Dintg_=$(intg_) -Dnev_=$(nev_) -Dme_=$(me_) \
			  -DV0_=$(V0_) -DB_i_=$(B_i_) -DB_f_=$(B_f_) -Dnum_puntos_B_=$(num_puntos_B_)

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
