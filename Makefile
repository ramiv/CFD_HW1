FC=gfortran
#FCFLAGS= -g -Wall -fdefault-real-8 -ffree-line-length-0 -fbacktrace
FCFLAGS= -g -fdefault-real-8 -ffree-line-length-0 -fbacktrace
OUTNAME = HW1.exe

SRCFILES_1 = $(wildcard src/*.f90)
SRCFILES = $(SRCFILES_1:src/%=%)

#obj files
OBJFILES = $(SRCFILES:.f90=.o)

OBJDIR = obj
SRCDIR = src

.PHONY: all

all: $(OBJDIR) $(OUTNAME) TAGS

TAGS: 
	ctags -R -f .

$(OBJDIR):
	mkdir $(OBJDIR)

$(OUTNAME): $(OBJFILES)
	@echo obj:$(OBJFILES)
	$(FC) $(FCFLAGS) -o $(OUTNAME) -I$(OBJDIR) -J$(OBJDIR) $(addprefix $(OBJDIR)/,$(OBJFILES))

HW1.o: HW1.f90 inmod.o mathmod.o

mathmod.o: mathmod.f90 inmod.o outmod.o diagal.o

inmod.o: inmod.f90

outmod.o: outmod.f90 inmod.o

%.mod: %.f90
	@echo compiling obj
	@echo target= $@
	@echo srcs = $<
	$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -o $(OBJDIR)/$($@:mod=.o) -c $<
%.o: %.f90
	$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -o $(OBJDIR)/$@ -c $<

%.f90:
	@echo compiling src
	@echo target= $@
	@echo output $(addprefix $(OBJDIR)/,$($@:.f90=.o))
	$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -o $(addprefix $(OBJDIR)/,$($@:.f90=.o)) -c $@

vpath %.f90 $(SRCDIR)
vpath %.o $(OBJDIR)
vpath %.mod $(OBJDIR)
