CC = gcc
CCFLAGS = -Wall -O2
LIBS = -lm -lgsl -lgslcblas

CC_COMPILE = $(CC) $(CCFLAGS) -c -fPIC
CC_LOAD = $(CC) $(CCFLAGS)

.c.o:
	$(CC_COMPILE) $*.c

EXE = gyst
all: $(EXE)

SRCS = main.c read_data.c sheets.c malloc.c eigen.c metric.c characterization.c

OBJS = main.o read_data.o sheets.o malloc.o eigen.o metric.o characterization.o

INCS = declarations.h definitions.h

$(OBJS) : $(INCS) makefile

$(EXE): $(OBJS) $(INCS) makefile
	$(CC_LOAD) $(OBJS) $(LIBS) -o $(EXE)

clean:
	/bin/rm -f *.o *_s $(EXE)

newrun:
	/bin/rm -rf run
	/bin/mkdir run
	/bin/mkdir run/images run/sheets_disc run/sheets_outflows run/sheets_jet \
	run/jcs_files_disc run/jcs_files_outflows run/jcs_files_jet
