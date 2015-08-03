#FC = mpif90 -mcmodel=medium -i-dynamic
# -check bound -check uninit
# FC = mpif90
FC = mpif90

# FFLAGS = -O2
# FFLAGS = -O2 -axW -132

F1 = compton2d.f \
	reader.f \
	setup2d.f \
	xec2d.f \
	imcgen2d.f \
	volume2d.f \
	gamma1_2d.f \
	nontherm2d.f \
	imcsurf2d_para.f \
	planck2d.f \
	pp2d.f \
	imctrk2d.f \
	census2d.f \
	compb_2d.f \
	comtot2d.f \
	imcdate2d.f \
	ref_matrix.f \
	imcleak2d.f \
	graphics2d.f \
	imcvol2d_para.f \
	imcfield2d.f \
        imcredist.f \
	icloss2d.f \
	update2d.f \
	rand.f \
	fp_mpi.f \
	surf_mpi.f \
	vol_mpi.f \
        write_record.f \
        read_record.f

O1 = compton2d.o \
	reader.o \
	setup2d.o \
	xec2d.o \
	imcgen2d.o \
	volume2d.o \
	gamma1_2d.o \
	nontherm2d.o \
	imcsurf2d_para.o \
	planck2d.o \
	pp2d.o \
	imctrk2d.o \
	census2d.o \
	compb_2d.o \
	comtot2d.o \
	imcdate2d.o \
	ref_matrix.o \
	imcleak2d.o \
	graphics2d.o \
	imcvol2d_para.o \
	imcfield2d.o \
        imcredist.o \
	icloss2d.o \
	update2d.o \
	rand.o \
	fp_mpi.o \
	surf_mpi.o \
	vol_mpi.o \
        write_record.o \
        read_record.o

compton2d : $(F1) $(O1)
	$(FC) -o compton2d $(O1)

.SUFFIXES: .f .o

.f.o:
	$(FC) -c -o $*.o $<

clean: 
	rm -f $(O1)

