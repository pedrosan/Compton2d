      subroutine IC_loss
      implicit none
      include 'general.pa'
      include 'commonblock.f'
c
c
c     Inverse-Compton loss of a single electron,
c       using full Klein-Nishina cross section,
c               no relativistic limits
c
c
c
      integer i_el, i_ph
c
      double precision gamma, epsilon
      double precision beta, F, A, z1, z2
      double precision f1, f2
c 
c
c        A = (c*pi*r_0)/2
c
      A = 3.7419d-15
c
      do 800 i_el = 1, num_nt
         gamma = 1.d0 + gnt(i_el)
         beta = dsqrt(1.d0 - 1.d0/(gamma*gamma))
c
         do 700 i_ph = 1, nphfield
            epsilon = 1.957d-3*E_field(i_ph)
c
            if ((gamma*epsilon).lt.1.d-2) then
               F_IC(i_el, i_ph) = 2.66d-14*epsilon
     1                           *(gamma*gamma - 1.d0)
            else
               z1 = epsilon*gamma*(1.d0 + beta)
               z2 = epsilon/(gamma*(1.d0 + beta))
c
               F = gamma*(f1(z1) - f1(z2)) 
     1           - epsilon*(f2(z1) - f2(z2))
c
               F_IC(i_el, i_ph) = A*F/(((epsilon*gamma)**2)*beta)
            endif
c
 700     continue
 800  continue
c
      open (unit=22, file='icloss.dat', status='unknown')
 850  format (e14.7,1x,e14.7,1x,e14.7,1x,e14.7)
c
      do 880 i_el =1, num_nt
        do 870 i_ph = 1, nphfield
           write(22, 850) gnt(i_el), E_field(i_ph), F_IC(i_el, i_ph)
 870    continue ! Xuhui Chen 11/7/2008
 880  continue

c      do 870 i_el = 1, num_nt
c         write(22, 850) gnt(i_el), F_IC(i_el, 1), F_IC(i_el, 50),
c     1                 F_IC(i_el, 100) 
c 870  continue
c
      close(22)
c
      return
      end
c
c
c
      double precision function f1(z)
      double precision z
c
      double precision f_Li, sd1, sd2, sum, y
c
      y = 1.d0 + 2.d0*z
      sd1 = (z + 6.d0 + 3.d0/z)*dlog(y)
      sd2 = ((2.2d1/3.d0)*(z**3) + 24.d0*(z**2)
     1    + 18.d0*z + 4.d0)/(y**2)
      sum = sd1 - sd2 - 2.d0 + 2.d0*f_Li(z)
      f1 = sum
c
      return
      end
c
c
c
      double precision function f2(z)
      implicit none
      double precision z
c
      double precision f_Li, sd1, sd2, sum, y
c
      y = 1.d0 + 2.d0*z
      sd1 = (z + (3.1d1/6.d0) + 5.d0/z + 1.5d0/(z**2))*dlog(y)
      sd2 = ((2.2d1/3.d0)*(z**3) + 28.d0*(z**2)
     1    + (1.03d2/3.d0)*z + 1.7d1 + 3.d0/z)/(y**2)
      sum = sd1 - sd2 - 2.d0 + f_Li(z)
      f2 = sum
c
      return
      end
c
c
c
c
      double precision function f_Li(z)
      implicit none
      double precision z
c
      double precision sum, sd, nn, y
c
      sum = 0.d0
      nn = 1.d0
      y = 1.d0 + 2.d0*z
      sum = dlog(y)*(5.d-1*dlog(y) - dlog(2.d0*z))
 200  if (nn.lt.1.5) then
         sd = 1.d0/y
      else
         sd = sd*(((nn - 1.d0)/nn)**2)/y
      endif
      nn = nn + 1.d0
      sum = sum + sd
      if (dabs(sd/sum).gt.1.d-10) goto 200
c
      f_Li = sum
      return
      end
c
