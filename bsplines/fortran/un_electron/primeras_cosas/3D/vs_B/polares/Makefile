# FILE
BIN		= cl1e

# MODULES
MOD1	= mod_prec
MOD2	= mod_cuad
MOD3	= mod_matrices
MOD4	= mod_eig

# COMPILER
CXX		= gfortran

# FLAGS
FLAGS			= -cpp -O3 -Wall -Wextra -fcheck=all
LDFLAGS		= -llapack

# VALUES
ZMIN_						= -1000.0
ZMAX_						= 1000.0
TIPO_						= 2
BETA_						= 0.0065d0
KORD_						= 5
L_INTERVAL_			= 15
N_CUAD_					= 200
NEV_						= 60
ME_							= 0.063d0
SIGMA_					= 15.0d0
V0_							= 0.00d0
B_CAMPO_I_			= 0.0d0
B_CAMPO_F_			= 100.0d0
NUM_PUNTOS_B_		= 100


# PARAMETER
PARAMETER	= -DZMIN_=$(ZMIN_) -DZMAX_=$(ZMAX_) -DTIPO_=$(TIPO_) -DBETA_=$(BETA_) \
					  -DKORD_=$(KORD_) -DL_INTERVAL_=$(L_INTERVAL_) -DN_CUAD_=$(N_CUAD_) \
						-DNEV_=$(NEV_) -DME_=$(ME_) -DSIGMA_=$(SIGMA_) -DV0_=$(V0_) -DB_CAMPO_I_=$(B_CAMPO_I_) \
						-DB_CAMPO_F_=$(B_CAMPO_F_) -DNUM_PUNTOS_B_=$(NUM_PUNTOS_B_)

$(BIN): clean $(BIN).f90 mods
	$(CXX) $(FLAGS) -o $(BIN) $(BIN).f90 $(MOD1).o $(MOD2).o $(MOD3).o $(MOD4).o $(LDFLAGS) $(PARAMETER)

run: $(BIN)
	./$(BIN) # > 'salida.dat'

mods:
	$(CXX) -c $(MOD1).f90 $(MOD2).f90 $(MOD3).f90 $(MOD4).f90

clean:
	rm -f $(BIN) *.mod *.o
