# Binary file
BIN					= cl2e

# Flags
CFLAGS			= -O3 -g -cpp -std=gnu99 -Wall -Wextra -fopenmp
LDFLAGS			= -llapack -lm

# Default Values
RMIN				= -800.0
RMAX				= 800.0
L						= 60
KORD				= 5
BETA				= 0.0065
TYPE				= 2
INTG				= 60
NEV					= 40
ME					= 0.063
V_POZO			= 0.05
SIGMA				= 20.0
B_CAMPO			= 80.0
LAMBDA_I		= 0.0
LAMBDA_F		= 1.0
NUM_LAMBDA	= 20

# Simulation Parameters
PARAMETERS	= -DRMIN=$(RMIN) -DRMAX=$(RMAX) -DL=$(L) -DKORD=$(KORD) -DBETA=$(BETA) -DTYPE=$(TYPE) -DINTG=$(INTG) -DNEV=$(NEV) \
							-DME=$(ME) -DV_POZO=$(V_POZO) -DSIGMA=$(SIGMA) -DB_CAMPO=$(B_CAMPO) \
							-DLAMBDA_I=$(LAMBDA_I) -DLAMBDA_F=$(LAMBDA_F) -DNUM_LAMBDA=$(NUM_LAMBDA)
# Compilers
CC					= gcc
LINKER			= gcc

# Files
C_SOURCES		= $(BIN).c
HEADERS			=
C_OBJS			= $(patsubst %.c, %.o, $(C_SOURCES))

# Rules

$(BIN): clean $(C_OBJS) $(HEADERS)
	$(LINKER) $(CFLAGS) -o $(BIN) $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CC) -c $(C_SOURCES) $(CFLAGS) $(INCLUDES) $(PARAMETERS)

run: $(MP) $(BIN)
	./$(BIN) #> $(OFILE) &

clean:
	rm -f $(BIN) *.o *.mod
