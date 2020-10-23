FC        = gfortran
RM        = rm -rf
LINKER    = $(FC)
LIBS      = 

ifeq ($(MAKECMDGOALS),debug) 
  COPT    = -c -O0 -cpp -fno-range-check -g -Wall -fbounds-check 
  CFLAG   =
  LOPT    =
else
  COPT    = -c -O3 -cpp -fno-range-check 
  CFLAG   =
  LOPT    = 
endif

EXECDIR   = ./

VPATH     = src

APNAME    = 4dvar 

TARGET    = mod_variables.o sub_check_foa.o m_random3.o sub_4dvar.o         \
            sub_measure.o sub_plot3d.o sub_bc_outer.o sub_bc_wall.o         \
            sub_hsmac_fwd.o sub_hsmac_tlm.o sub_hsmac_adj.o sub_initial.o   \
            sub_rhs3rd_fwd.o sub_rhs3rd_tlm.o sub_rhs3rd_adj.o sub_utils.o  \
            sub_bgfilter.o NS2D_4DVar.o 

CLNTARGET = *.o *.mod *.d *.pc *.L $(APNAME)


$(APNAME) : $(TARGET)
	$(LINKER) $(CFLAG) -o $(EXECDIR)$(APNAME) $(TARGET) $(LIBS)

%.o : %.F90
	$(FC) $(COPT) -c -o $@ $<
%.o : %.f90
	$(FC) $(COPT) -c -o $@ $<
%.o : %.f
	$(FC) $(COPT) -c -o $@ $<

debug : $(TARGET)
	$(LINKER) $(CFLAG) -o $(EXECDIR)$(APNAME) $(TARGET) $(LIBS)

clean :
	$(RM) $(CLNTARGET)

realclean :
	$(RM) $(CLNTARGET) $(APNAME) *~
