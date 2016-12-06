# Binary file
BIN			= 1eclvsB

# Flags
CFLAGS			= -O3 -cpp -Wall -Wextra
LDFLAGS			= -llapack #-mkl=sequential -lmkl_lapack95_lp64 -lpthread -lm

# Default Values
zmin_						= -1000.d0
zmax_						= 1000.d0
lum_						= 0
c_							= 0.d0
gamma_					= 0.065d0
kord_						= 7
l_interval_			= 150
intg_						= 100
nev_						= 50
V1_							= 0.37d0
V2_ 						= 0.108844d0
r0_ 						= 7.0d0
az_ 						= 7.0d0
bz_ 						= 2.5d0
B_campo_i_			= 1.d0
B_campo_f_			= 50.d0
num_puntos_B_		= 50
me_							= 0.041d0


# Simulation Parameters
PARAMETERS		= -Dzmin_=$(zmin_) -Dzmax_=$(zmax_) -Dlum_=$(lum_) -Dc_=$(c_) \
			  -Dgamma_=$(gamma_) -Dkord_=$(kord_) -Dl_interval_=$(l_interval_) -Dintg_=$(intg_) -Dnev_=$(nev_) -DV1_=$(V1_) -DV2_=$(V2_) \
			  -Dr0_=$(r0_) -Daz_=$(az_) -Dbz_=$(bz_) -DB_campo_i_=$(B_campo_i_) -DB_campo_f_=$(B_campo_f_) \
				-Dnum_puntos_B_=$(num_puntos_B_) -Dme_=$(me_)

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
