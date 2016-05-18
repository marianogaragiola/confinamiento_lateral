# FILE
BIN		= cl2e5

# COMPILER
CXX		= gfortran

# FLAGS
FLAGS		= -cpp -ggdb -O3
LDFLAGS		= -llapack

# VALUES
me_		= 0.063d0
V0_		= 0.05d0
B_campo_i_	= 5.d0
B_campo_f_	= 100.d0
N_B_		= 95
sigma_		= 20.d0
lambda_		= 0.69d0
omega_		= 0.030 #0213877d0

# PARAMETER
PARAMETER	= -Dme_=$(me_) -DV0_=$(V0_) -DB_campo_i_=$(B_campo_i_) -DB_campo_f_=$(B_campo_f_) \
		  -DN_B_=$(N_B_) -Dsigma_=$(sigma_) -Dlambda_=$(lambda_) -Domega_=$(omega_)

$(BIN): clean $(BIN).f90
	$(CXX) $(FLAGS) -o $(BIN) $(BIN).f90 $(LDFLAGS) $(PARAMETER)

run: $(BIN)
	./$(BIN) # > 'salida.dat'

clean:
	rm -f $(BIN) *.mod
