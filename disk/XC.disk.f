      program main
      implicit none
      
      double precision h, c, k, sigma, pi, nu_alpha
      parameter (sigma=5.67d-5)
      parameter (pi=3.14159d0)
      parameter (h=6.62618d-27)
      parameter (c=2.99792d+10)
      parameter (k=1.38d-16)
      parameter (nu_alpha=2.47d15)
      integer nphmax, nph_tave, i, j
      parameter (nphmax=500) 
      parameter (nph_tave=5402)
      double precision nu, dnu, T_disk, T_blr, T_ir, E, Eph(nphmax),
     1 Fph_disk(nphmax), Fph_blr(nphmax), Fph_tave(nphmax),
     1 Fph_ir(nphmax), E_min, E_max, pnth, Gam, bet, E_tave(nph_tave)
      double precision Utherm_disk, Utherm_blr, Unth_disk, Unth_blr,
     1 Utherm_ir, Unth_ir, Umax, y, L_disk
      double precision Fph_beam(nphmax), Fph_unbeam(nph_tave)
      
      Gam=15.d0
      bet=sqrt(1.d0-1.d0/Gam**2)
      T_disk = 3.d4
      T_blr = 1.5*Gam*nu_alpha*h/3.93/k !from GG09
      T_ir =367.d0*Gam !from GG09
      E_min = 5d7 !5.d-4
      E_max = 5d8 !1.d-1
      pnth = 1.d0
      E = 1.d-7 ! photon energy in keV
      nu = E*1.602d-9/h
      dnu = 10**(log10(1d10/E)/nphmax)
      E=E*sqrt(dnu) ! E is the median energy of nphmax energy bins
      Umax =0.d0
      do 100 i=1, nphmax
        if(E.le.E_min)then
          Utherm_disk = 9.d62* 2*h*nu**3/c**2/(exp(h*nu/k/T_disk)-1.d0)
     1             /(sigma/pi*T_disk**4)
          Utherm_blr = 1.d44* 2*h*nu**3/c**2/(exp(h*nu/k/T_blr)-1.d0)
     1             /(sigma/pi*T_blr**4)
          Utherm_ir = 1.d44* 2*h*nu**3/c**2/(exp(h*nu/k/T_ir)-1.d0)
     1             /(sigma/pi*T_ir**4)

c          if(Utherm.gt.Umax)Umax = Utherm
c          Eph(i) = E
          Fph_disk(i) = Utherm_disk
          Fph_blr(i) = Utherm_blr
          Fph_ir(i) = Utherm_ir
        else
          y = E/E_max
          if(y.lt.1.d2)then
            Unth_disk = Utherm_disk*(E/E_min)**(-pnth)/exp(y)
            Unth_blr = Utherm_blr*(E/E_min)**(-pnth)/exp(y)
            Unth_ir = Unth_ir*(E/E_min)**(-pnth)/exp(y)
          else
            Unth_disk = 0.d0
            Unth_blr = 0.d0
            Unth_ir = 0.d0
          endif
          Fph_disk(i) = Unth_disk
          Fph_blr(i) = Unth_blr
          Fph_ir(i) = Unth_ir
        endif
        Eph(i) = E
        nu = nu*dnu
        E = E*dnu
 100  continue
      open(unit=6, file='tave_ub.dat', status='unknown')
      open(unit=5, file='tave_f.dat', status='unknown')
      open(unit=4, file='blackbody.in', status='unknown')
      do i=1,nph_tave
        read(6,fmt='(e14.7, 1x, e14.7)')E_tave(i), Fph_unbeam(i)
      enddo
      E_tave(:)=10**E_tave(:)*h/1.602d-9
      Fph_unbeam(:)=10**Fph_unbeam(:)
c     eq(4) in Tavecchio et al. 2008, but changed all nu into E. 
c     And F=U*c under the assumption that the beamed radiation all comes 
c     from theta=0
      do i=1,nphmax
        Fph_beam(i)=0
        do j=1,nph_tave-1
          if(E_tave(j).gt.Eph(i)/Gam/(1.d0+bet).and.
     1     E_tave(j).le.Eph(i)/Gam)then
           Fph_beam(i)=Fph_beam(i)+Fph_unbeam(j)/E_tave(j)**3*
     1     (E_tave(j+1)-E_tave(j))
          endif
        enddo
        Fph_beam(i)=Fph_beam(i)*2*pi*Eph(i)**2/Gam/bet
      enddo
      do 200 i=1,nphmax
       read(5,fmt= '(e14.7, 1x, e14.7)')E, Fph_tave(i)
       write(4,fmt= '(e14.6,3(1x,e14.6))') Eph(i),
     1 dmax1(1.d-30, Fph_disk(i)), dmax1(1.d-30,Fph_beam(i)),
     1 dmax1(1.d-30, Fph_ir(i))
 200  continue
      L_disk = 0
      do i=1,nphmax-1
         L_disk = L_disk+Fph_disk(i)*(Eph(i+1)-Eph(i))
      enddo
      write(*,*)'L_disk=',L_disk
      close(4)
      close(5)
      end


