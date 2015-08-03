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

      double precision Gamma, beta
      double precision nu, dnu, E, E_min, E_max, y, pnth
      double precision T_disk, T_blr, T_ir
      double precision Utherm_disk, Utherm_blr, Unth_disk, Unth_blr, Utherm_ir, Unth_ir, Umax, L_disk

      double precision Eph(nphmax), Fph_disk(nphmax), Fph_blr(nphmax), Fph_ir(nphmax)
      double precision Fph_beamed_BLR(nphmax), Fph_tave(nphmax)
      double precision E_tave(nph_tave), Fph_unbeamed(nph_tave)  


      write(6,100) 
      read(5,*) Gamma
100   format(' Lorentz factor : ',$)

      ! Gamma= 15.d0
      beta = sqrt(1.d0 - 1.d0/Gamma**2)

      T_disk = 3.d4
      T_blr  = 1.5*Gamma*nu_alpha*h/3.93/k !from GG09
      T_ir   = 367.d0*Gamma                !from GG09

      E_min = 5d7 !5.d-4
      E_max = 5d8 !1.d-1
      pnth = 1.d0
      E = 1.d-7 ! photon energy in keV
      nu = E*1.602d-9/h
      dnu = 10**(log10(1d10/E)/nphmax)
      E = E*sqrt(dnu) ! E is the median energy of nphmax energy bins

      Umax =0.d0
      do i=1, nphmax

         if(E .le. E_min)then
            Utherm_disk = 9.0d62*2*h*nu**3/c**2/(exp(h*nu/k/T_disk)-1.d0)/(sigma/pi*T_disk**4)
            Utherm_blr  = 1.0d44*2*h*nu**3/c**2/(exp(h*nu/k/T_blr) -1.d0)/(sigma/pi*T_blr**4)
            Utherm_ir   = 1.0d44*2*h*nu**3/c**2/(exp(h*nu/k/T_ir)  -1.d0)/(sigma/pi*T_ir**4)

!            if(Utherm.gt.Umax)Umax = Utherm
!            Eph(i) = E
            Fph_disk(i) = Utherm_disk
            Fph_blr(i)  = Utherm_blr
            Fph_ir(i)   = Utherm_ir
         else
            y = E/E_max
            if(y .lt. 1.0d2)then
              Unth_disk = Utherm_disk*(E/E_min)**(-pnth)/exp(y)
              Unth_blr  = Utherm_blr*(E/E_min)**(-pnth)/exp(y)
              Unth_ir   = Unth_ir*(E/E_min)**(-pnth)/exp(y)
            else
              Unth_disk = 0.d0
              Unth_blr  = 0.d0
              Unth_ir   = 0.d0
            end if
            Fph_disk(i) = Unth_disk
            Fph_blr(i)  = Unth_blr
            Fph_ir(i)   = Unth_ir
         end if
         Eph(i) = E
         nu     = nu*dnu
         E      = E*dnu

      end do

      !open(unit=26, file='tave_ub.dat', status='old')
      !open(unit=25, file='tave_f.dat', status='old')
      open(unit=26, file='tavecchio_Uext.dat', status='old')
      open(unit=25, file='tavecchio_Uext_beamed.dat', status='old')
      open(unit=24, file='blackbody.in', status='unknown')

      do i=1,nph_tave
         ! read(26,fmt='(e14.7, 1x, e14.7)') E_tave(i), Fph_unbeamed(i)
         read(26,fmt='(5x,f8.4,5x,f8.4)') E_tave(i), Fph_unbeamed(i)
      end do
      E_tave(:)       = 10**E_tave(:)*h/1.602d-9
      Fph_unbeamed(:) = 10**Fph_unbeamed(:)
      ! eq(4) in Tavecchio et al. 2008, but changed all nu into E. 
      ! And F=U*c under the assumption that the beamed radiation all comes 
      ! from theta=0, mu=1

      do i=1,nphmax
         Fph_beamed_BLR(i)=0
         do j=1,nph_tave-1
            if( E_tave(j) .gt. Eph(i)/Gamma/(1.d0+beta) .and. E_tave(j) .le. Eph(i)/Gamma) then
               Fph_beamed_BLR(i) = Fph_beamed_BLR(i) + Fph_unbeamed(j)/E_tave(j)**3.*(E_tave(j+1)-E_tave(j))
            end if
         end do
         Fph_beamed_BLR(i) = Fph_beamed_BLR(i)*2*pi*Eph(i)**2/Gamma/beta
      end do

      do i=1,nphmax
       read(25,fmt= '(1x,e13.7,2x,e13.7)') E, Fph_tave(i)
       write(24,fmt= '(e14.6,4(1x,e14.6))') Eph(i),
     &                                      dmax1(1.d-30,Fph_disk(i)), 
     &                                      dmax1(1.d-30,Fph_beamed_BLR(i)),
     &                                      dmax1(1.d-30,Fph_ir(i)),
     &                                      dmax1(1.d-30,Fph_blr(i))
      end do


      L_disk = 0
      do i=1,nphmax-1
         L_disk = L_disk + Fph_disk(i)*(Eph(i+1)-Eph(i))
      enddo
      write(6,*)' Integrated L_disk = ',L_disk
      close(24)
      close(25)
      close(26)

      end

