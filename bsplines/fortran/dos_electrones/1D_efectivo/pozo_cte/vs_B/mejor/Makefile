# FILE
BIN		= cl2e

# MODULES
MOD1	= mod_prec
MOD2	= mod_cuad
MOD3	= mod_matrices
MOD4	= mod_eig

# COMPILER
CXX		= gfortran

# FLAGS
FLAGS			= -cpp -O0 -Wall -Wextra -fcheck=all
LDFLAGS		= -llapack

# VALUES
ZMIN_						= -1000.0
ZMAX_						= 1000.0
TIPO_						= 2
BETA_						= 0.0065d0
KORD_						= 7
L_INTERVAL_			= 100
N_CUAD_					= 100
NEV_						= 50
ME_							= 0.041d0
AZ_							= 7.d0
BZ_							= 2.5d0
R0_							= 7.d0
V1_							= 0.37d0
V2_							= 0.108844d0
ETA_						= 1.d0
DELTA_					= 0.1d0
BCAMPO_I_				= 1.d0
BCAMPO_F_				= 50.d0
NUM_PUNTOS_B_		= 50

# PARAMETER
PARAMETER	= -DZMIN_=$(ZMIN_) -DZMAX_=$(ZMAX_) -DTIPO_=$(TIPO_) -DBETA_=$(BETA_) \
					  -DKORD_=$(KORD_) -DL_INTERVAL_=$(L_INTERVAL_) -DN_CUAD_=$(N_CUAD_) \
						-DNEV_=$(NEV_) -DME_=$(ME_) -DAZ_=$(AZ_) -DBZ_=$(BZ_) -DR0_=$(R0_) \
						-DV1_=$(V1_) -DV2_=$(V2_) -DETA_=$(ETA_) -DDELTA_=$(DELTA_) \
						-DBCAMPO_I_=$(BCAMPO_I_) -DBCAMPO_F_=$(BCAMPO_F_) \
						-DNUM_PUNTOS_B_=$(NUM_PUNTOS_B_)

$(BIN): clean $(BIN).f90 mods
	$(CXX) $(FLAGS) -o $(BIN) $(BIN).f90 $(MOD1).o $(MOD2).o $(MOD3).o $(MOD4).o $(LDFLAGS) $(PARAMETER)

run: $(BIN)
	./$(BIN) # > 'salida.dat'

mods:
	$(CXX) -c $(MOD1).f90 $(MOD2).f90 $(MOD3).f90 $(MOD4).f90

clean:
	rm -f $(BIN) *.mod *.o
