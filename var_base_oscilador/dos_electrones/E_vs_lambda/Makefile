# FILE
BIN		= cl2e4

# COMPILER
CXX		= gfortran

# FLAGS
FLAGS		= -cpp -ggdb -O3
LDFLAGS		= -llapack

# VALUES
me_		= 0.063d0
V0_		= 0.05d0
B_campo_	= 80.d0
omega_		= 0.00001d0
lambda_i_	= 0.d0
lambda_f_	= 1.d0
N_lambda_	= 10
sigma_		= 20.d0

# PARAMETER
PARAMETER	= -Dme_=$(me_) -DV0_=$(V0_) -DB_campo_=$(B_campo_) -Domega_=$(omega_) -Dlambda_i_=$(lambda_i_) -Dlambda_f_=$(lambda_f_) \
		  -DN_lambda_=$(N_lambda_) -Dsigma_=$(sigma_)

$(BIN): clean $(BIN).f90
	$(CXX) $(FLAGS) -o $(BIN) $(BIN).f90 $(LDFLAGS) $(PARAMETER)

run: $(BIN)
	./$(BIN) # > 'salida.dat'

clean:
	rm -f $(BIN) *.mod
