# dimer


CFLAGS = -c -Ofast -Wall
SRCPATH = src/
OBJPATH = obj/
DATAPATH = data/
FIGPATH = figs/
SOURCES = $(wildcard $(SRCPATH)*.c)
OBJECTS = $(patsubst $(SRCPATH)%.c, $(OBJPATH)%.o, $(SOURCES))
EXECUTABLE = $(OBJPATH)prog
RUN = mpirun -np # run command of programm 'prog'
N_PROC = 2 # number of the processors


all:
	$(RUN) $(N_PROC) $(EXECUTABLE)
	python3 plot.py


all: $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	mpicc $(OBJECTS) -lm -o $@


$(OBJPATH)%.o: $(SRCPATH)%.c
	mpicc $(CFLAGS) $< -o $@


$(OBJECTS): $(SRCPATH)header.h makefile


# command $make clean
clean:
	rm -rf $(OBJPATH)
	rm -rf $(DATAPATH)
	rm -rf MPI_data/
	rm -rf $(FIGPATH)
	mkdir $(OBJPATH)
	mkdir $(DATAPATH)
	mkdir MPI_data/
	mkdir $(FIGPATH)
