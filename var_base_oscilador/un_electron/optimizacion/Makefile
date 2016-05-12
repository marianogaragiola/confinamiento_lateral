# FILE
BIN		= cl1e2

# COMPILER
CXX		= gfortran

# FLAGS
FLAGS		= -cpp -ggdb -O3
LDFLAGS		= -llapack

# VALUES
me_		= 0.063d0
V0_		= 0.05d0
omega_i_	= 0.d0
omega_f_	= 1.d0
N_omega_	= 100
sigma_		= 20.d0
B_campo_	= 80.d0

# PARAMETER
PARAMETER	= -Dme_=$(me_) -DV0_=$(V0_) -Domega_i_=$(omega_i_) -Domega_f_=$(omega_f_) \
		  -DN_omega_=$(N_omega_) -Dsigma_=$(sigma_) -DB_campo_=$(B_campo_)

$(BIN): clean $(BIN).f90
	$(CXX) $(FLAGS) -o $(BIN) $(BIN).f90 $(LDFLAGS) $(PARAMETER)

run: $(BIN)
	./$(BIN) # > 'salida.dat'

clean:
	rm -f $(BIN) *.mod
