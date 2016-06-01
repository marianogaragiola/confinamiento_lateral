# Binary file
BIN					= cl1e

# Flags
CFLAGS			= -O3 -cpp -std=gnu99 -Wall -Wextra -fopenmp
LDFLAGS			= -llapack -lm

# Default Values
RMIN				= -200.0
RMAX				= 200.0
L						= 70
KORD				= 5
INTG				= 60
NEV					= 40
ME					= 0.063
V_POZO			= 0.05
SIGMA				= 20.0
B_CAMPO_I		= 10.0
B_CAMPO_F		= 100.0
NUM_B				= 90

# Simulation Parameters
PARAMETERS	= -DRMIN=$(RMIN) -DRMAX=$(RMAX) -DL=$(L) -DKORD=$(KORD) -DINTG=$(INTG) -DNEV=$(NEV) \
							-DME=$(ME) -DV_POZO=$(V_POZO) -DSIGMA=$(SIGMA) -DB_CAMPO_I=$(B_CAMPO_I) \
							-DB_CAMPO_F=$(B_CAMPO_F) -DNUM_B=$(NUM_B)
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
