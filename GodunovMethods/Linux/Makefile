BASE = $(shell /bin/pwd)
EXED = $(BASE)/$(EXECTBL)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src

vpath %.o $(OBJD)
vpath %.f90 $(SRCD)

EXECTBL  = conservative

SRCS =	conservative.f90

OBJS =	conservative.o

LIBS =	

CC = cc
CFLAGS = -O
FC = pgf77
FFLAGS = -O
# insert below the f90 compiler available, eg F90 = pgf90 or F90 = ifort
#F90 = g95
F90 = /opt/intel/bin/ifort
#
F90FLAGS = -O
#LDFLAGS = -s
LDFLAGS = 

all: $(EXECTBL)

$(EXECTBL): $(OBJS)
	@echo '' 
	@echo '============ ...building the executable ============'
	@echo ''
	cd $(OBJD); $(F90) $(F90FLAGS) -o ../$(EXECTBL) $(OBJS)
	@echo ''
	@echo '============ ...Done! ============'

%.o $(OBJD)/%.o: %.f90
	@echo '' 
	@echo '============ ...building' $*.o '============'
	@echo ''
	cd $(OBJD); $(F90) -c $(F90FLAGS) $(SRCD)/$(*).f90

flush:
	rm -f $(EXECTBL) *.dat $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc*

clean:
	rm -f $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc* $(SRCD)/*.o $(EXECTBL)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

# DO NOT DELETE
