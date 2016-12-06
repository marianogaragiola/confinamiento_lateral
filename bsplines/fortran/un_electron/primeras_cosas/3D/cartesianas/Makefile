# Binary file
BIN					= cl1e

# Flags
CFLAGS			= -O3 -Wall -Wextra -cpp
LDFLAGS			= -llapack

# Default Values
zmin_						= -1000.d0
zmax_						= 1000.d0
z0_							= 0.d0
sigma_					= 15.d0
lum_						= 0
c_							= 0.d0
beta_						= 0.0065d0
kord_						= 5
l_							= 10
intg_						= 100
nev_						= 30
me_							= 0.063d0
V0_i_						= 0.d0
V0_f_						= 0.03d0
num_puntos_V0_	= 30

# Simulation Parameters
PARAMETERS		= -Dzmin_=$(zmin_) -Dzmax_=$(zmax_) -Dz0_=$(z0_) -Dsigma_=$(sigma_) -Dlum_=$(lum_) -Dc_=$(c_) \
							  -Dbeta_=$(beta_) -Dkord_=$(kord_) -Dl_=$(l_) -Dintg_=$(intg_) -Dnev_=$(nev_) -Dme_=$(me_) \
			 				  -DV0_i_=$(V0_i_) -DV0_f_=$(V0_f_) -Dnum_puntos_V0_=$(num_puntos_V0_)

# Compilers
CC					= gfortran
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
