SRC_Folder = src
MOD_Folder = mod


FC = gfortran
FCFLAGS = -g -Wall -fdefault-real-8 -ffree-line-length-0 -fbacktrace
MODFLAG = -J
INCFLAG = -I



target_try = try
target_all = HW1


help:
	@echo "$(SRC_Folder)"
	@echo "$(MOD_Folder)"
	@echo "$(MOD_Folder)/%.o:"

all: $(target_all)


$(target_try):	$(MOD_Folder)/tryMods.o
	$(FC) -o $@  $< $(MOD_Folder)/inmod.o $(MOD_Folder)/outmod.o $(MOD_Folder)/mathmod.o $(MODFLAG)$(MOD_Folder) 
	ctags -R .

$(target_all):	$(MOD_Folder)/HW1.o
	$(FC) -o $@  $< $(MOD_Folder)/inmod.o $(MOD_Folder)/outmod.o $(MOD_Folder)/mathmod.o $(MOD_Folder)/diagal.o $(INCFLAG)$(MOD_Folder) $(MODFLAG)$(MOD_Folder) 
	ctags -R .

# Default Rule
$(MOD_Folder)/%.o: $(SRC_Folder)/%.f90
	$(FC) -c $(FCFLAGS) $(INCFLAG)$(MOD_Folder) $(MODFLAG)$(MOD_Folder) $< -o $@

$(MOD_Folder)/tryMods.o: $(MOD_Folder)/inMod.o $(MOD_Folder)/outMod.o $(MOD_Folder)/mathmod.o

$(MOD_Folder)/HW1.o: $(MOD_Folder)/inMod.o $(MOD_Folder)/mathmod.o

$(MOD_Folder)/outmod.o: $(MOD_Folder)/inMod.o 
$(MOD_Folder)/mathmod.o: $(MOD_Folder)/inMod.o $(MOD_Folder)/diagal.o $(MOD_Folder)/outMod.o 

.PHONY:	clean

clean:
	rm -f $(target_try).exe $(MOD_Folder)/*.o $(MOD_Folder)/*.mod
