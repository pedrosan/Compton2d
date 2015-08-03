      implicit none
      include 'general.pa'
      include 'commonblock.f'
      character *30 fname
c
      integer,parameter :: ndisk = 10
       
c
      integer i, n_input, j, ijump
      double precision mu_disk(ndisk), dmu(ndisk)
      double precision r_d(ndisk),
     1  Ltot_disk, Ftot_blr, Ftot_ir, Ftot_blr_norm, Ftot_ir_norm
      double precision Id_file(nfmax,ndisk), F_como(nfmax,ndisk),
     1  L_disk(nfmax), F_blr(nfmax), F_ir(nfmax)
c
      double precision Isum
      double precision beta, dopp(ndisk)
c
c

c
      n_input = 25

      open(n_input, file=fname, status='unknown') ! Xuhui Chen 05/08
      i = 0
 100  i = i + 1
      read(unit=n_input, err=120, fmt=*) 
     1 E_file(i), L_disk(i), F_blr(i), F_ir(i)
c      L_disk(i)=L_disk(i)*max(1.d0,7.3d0*exp(-abs(time-8.d6)/2.d6)) ! Xuhui 1/13/10
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculation of radiative intensity from the accretion disk.
c     followed the calculation in Ghisellini & Tavecchio 2009.
c
c     d_jet is the distance from the jet to the accretion disk.
c     dopp is the doppler factor of the photons travel from the
c     accretion disk to the jet.
      do j=1,ndisk
        r_d(j)=R_disk*(j-0.5)/ndisk
      enddo
      beta = sqrt(1.d0-1.d0/g_bulk**2)
      mu_disk(:) = d_jet/sqrt(d_jet**2+r_d(:)**2)
      do j=1,ndisk
         dmu(j) = d_jet/sqrt(d_jet**2+(R_disk*(j-1)/ndisk)**2)-
     1   d_jet/sqrt(d_jet**2+(R_disk*j/ndisk)**2)
      enddo
      dopp(:) = g_bulk*(1.d0-beta*mu_disk(:))
      Id_file(i,:)=L_disk(i)/(4*pi**2*R_disk**2*mu_disk(:))
      F_como(i,:)=2*pi*dopp(:)*Id_file(i,:)*dmu(:)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((E_file(i).gt.0.d0).and.(i.le.(nfmax-1))) goto 100
 120  close(n_input)
      nfile = i - 1

c     The following procedure has assumed E_file to be evenly spaced
c     in logarithm. It moves the new spectrum to energy grids of its own
c      F_file(:)=1.d-29
c      do j=1,ndisk
c       if(dopp(j).ge.1.d0)then
c         do ijump=0,nfmax-1
c            if(E_file(ijump+1).ge.E_file(1)*dopp(j))exit
c         enddo
c         do i=1,nfmax-ijump
c            F_file(i+ijump)=F_file(i+ijump)+F_como(i,j)
c         enddo
c       else
c         do ijump=0,nfmax-1
c            if(E_file(nfmax-ijump).le.E_file(nfmax)*dopp(j))exit
c         enddo
c         do i=1,nfmax-ijump
c            F_file(i)=F_file(i)+F_como(i+ijump,j)
c         enddo
c       endif
c      enddo

c     Ltot_disk is the total power of the disk radiation.
c     E_file is the median of every energy grid, so sqrt() is divided for
c     the actual power.
      Ltot_disk = 0.d0
      Ftot_blr = 0.d0
      Ftot_ir = 0.d0
      do i=1,nfmax-1
        Ltot_disk = Ltot_disk+L_disk(i)*(E_file(i+1)-E_file(i))
        Ftot_blr = Ftot_blr+F_blr(i)*(E_file(i+1)-E_file(i))
        Ftot_ir = Ftot_ir+F_ir(i)*(E_file(i+1)-E_file(i))
      enddo
      Ltot_disk = Ltot_disk/sqrt(E_file(2)/E_file(1))
c     Ftot_blr is the total power of the BLR radiation of the input file.
c     Ftot_ir is the total power of the torus radiation of the input file.
      Ftot_blr = Ftot_blr/sqrt(E_file(2)/E_file(1))
      Ftot_ir = Ftot_ir/sqrt(E_file(2)/E_file(1))
c     Ftot_blr_norm and Ftot_ir_norm are the total power the BLR and torus 
c     radiation should have according
c     to Ghisellini & Madau 1996 (timed c for flux)
      Ftot_blr_norm = 17.d0/48.d0/pi*g_bulk**2*fr_blr*Ltot_disk/R_blr**2
      Ftot_ir_norm =  1.d0/4.d0/pi*g_bulk**2*fr_ir*Ltot_disk/R_ir**2
c     F_blr in the end is the normalized BLR and torus flux.
      F_blr(:) = F_blr(:)/Ftot_blr*Ftot_blr_norm
      F_ir(:) = F_ir(:)/Ftot_ir*Ftot_ir_norm


c************** Here it decides which of the above calculation counts*******
c      F_file(:)=F_file(:)+F_blr(:)
      F_file(:)=F_blr(:)

      if (nfile.lt.2) then
         write(*,*) 'Error in input spectrum file:'
         write(*,*) 'Less than 2 lines of input read!'
         write(*,*)
         stop
      endif
c
c      do 130 i = 1, nfile
c         write(*,*) i, E_file(i), F_file(i)
c 130  continue
c
      Isum = 0.d0
      do 150 i = 1,nfile-1
         alpha(i) = log(F_file(i+1)/F_file(i))
     1             /log(E_file(i+1)/E_file(i))
         a1(i) = alpha(i) + 1.d0
         if (a1(i).gt.2.d1) a1(i) = 2.d1
         if (a1(i).lt.-2.d1) a1(i) = -2.d1
         if (dabs(a1(i)).lt.1.d-3) then
            I_file(i) = F_file(i)*E_file(i)
     1                 *log(E_file(i+1)/E_file(i))
         else
            I_file(i) = F_file(i)*E_file(i)
     1                *((E_file(i+1)/E_file(i))**a1(i)
     2                  - 1.d0)/a1(i)
         endif
         Isum = Isum + I_file(i)
         P_file(i) = Isum
 150  continue
c
      do 200 i = 1, nfile-1
         P_file(i) = P_file(i)/Isum
 200  continue

      int_file0 = Isum
c
      return
      end

