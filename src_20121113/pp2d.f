c
c      These subroutines calculate the
c     differential pair production rate
c    
c
      subroutine pairprod(j, k)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer j, k
c
c
      integer i_p1, i_p2, i_g
c
      double precision gamma, sum, E1, eps1, E2, eps2, sum1, sd
      double precision f_pprod, emin
c
c
      if (j.eq.1.and.k.eq.1) 
     1   open(18, file='dn_e1.dat', status='unknown')
      if (j.eq.2.and.k.eq.1) 
     1   open(19, file='dn_e2.dat', status='unknown')
c
      do 500 i_g = 1, num_nt
         gamma = gnt(i_g) + 1.0d0
         sum = 0.d0
c
         do 400 i_p1 = 1, n_gg-1
            if (n_ph(i_p1, j, k).lt.1.d-5) goto 400
            E1 = E_gg(i_p1)
            eps1 = 1.957d-3*E1
            sum1 = 0.d0
            emin = dmax1((1.d0/eps1), (gamma + 1.d0 - eps1))
c
            do 300 i_p2 = 1, n_gg-1
               E2 = E_gg(i_p2)
               eps2 = 1.957d-3*E2
               if (eps2.lt.emin) goto 300
               if (n_ph(i_p2, j, k).lt.1.d-5) goto 300
c
               sd = f_pprod(eps1, eps2, gamma)
               sum1 = sum1 + sd*(E_gg(i_p2+1) - E2)
     1                      *n_ph(i_p2, j, k)/(eps2**2)
  300       continue
            sum = sum + sum1*(E_gg(i_p1+1) - E1)*n_ph(i_p1, j, k)
     1                 /(eps1**2)
  400    continue
c
         dn_pp(j, k, i_g) = 1.496d-14*sum
c
         if (j.eq.1.and.k.eq.1) then
           write(18, 600) gamma, dmax1(1.d-20, dn_pp(i_g, j, k))
         else if (j.eq.2.and.k.eq.1) then
           write(19, 600) gamma, dmax1(1.d-20, dn_pp(i_g, j, k))
         endif
c
  500 continue
c
  600 format (e10.3,1x,e10.3)
      if (j.eq.1.and.k.eq.1) then
         close(18)
      else if (j.eq.2.and.k.eq.1) then
         close(19)
      endif
c
      return
      end
c
c
c
      double precision function f_pprod(eps1, eps2, gamma)
      implicit none
      double precision eps1, eps2, gamma
c
      double precision f_inner, ecm_U, ecm_L, x
      double precision E, edagger, estar, Det
      double precision Det2, estar2, edagger2
c
      E = eps1 + eps2
      x = gamma*(E - gamma)
      Det2 = (x + 1.d0)**2 - E**2
      if (Det2.lt.0.d0) then
         f_pprod = 0.d0
         goto 900
      endif
      Det = dsqrt(Det2)
c
      estar2 = 5.d-1*(x + 1.d0 + Det)
      edagger2 = 5.d-1*(x + 1.d0 - Det)
      if ((estar2.lt.0.d0).or.(edagger2.lt.0.d0)) then
         f_pprod = 0.d0
         goto 900
      endif
      estar = dsqrt(estar2)
      edagger = dsqrt(edagger2)
c        
      ecm_U = dmin1(dsqrt(eps1*eps2), estar)
      ecm_L = dmax1(1.d0, edagger)
c     
      f_pprod = f_inner(ecm_U, eps1, eps2, gamma)
     1        - f_inner(ecm_L, eps1, eps2, gamma)
c
 900  return
      end
c
c
c
      double precision function f_inner(ecm, eps1, eps2, gamma)
      implicit none
      double precision ecm, eps1, eps2, gamma
c
      double precision E, f1, f12, H
c
      E = eps1 + eps2
c
      f12 = E**2 - 4.d0*(ecm**2)
      if (f12.lt.0.d0) then
        f_inner = 0.d0
        goto 900
      endif
      f1 = 2.5d-1*dsqrt(f12)
c
      f_inner = f1 + H(ecm, eps1, eps2, gamma)
     1             + H(ecm, eps2, eps1, gamma)
c
 900  return
      end
c
c
c
      double precision function H(ecm, eps1, eps2, gamma)
      implicit none
      double precision ecm, eps1, eps2, gamma
c
      double precision c, d, d2, I_pm, ee
c
      ee = eps1*eps2
      c = (eps1 - gamma)**2 - 1.d0
      d = eps1**2 + ee + gamma*(eps2 - eps1)
      d2 = ee + c*(ecm**2)
      if (d2.le.0.d0) then
         H = 0.d0
         goto 900
      endif
c
      if (dabs(c).gt.1.d-10) then
         H = -1.25d-1*ecm*(d/ee + 2.d0/c)/dsqrt(d2)
     1     + 2.5d-1*(2.d0 - (ee - 1.d0)/c)
     2             *I_pm(ecm, eps1, eps2, c)
     3     + 2.5d-1*dsqrt(d2)*(ecm/c + 1.d0/(ecm*ee))
      else
         H = ((ecm**3)/1.2d1 - 1.25d-1*ecm*d)/(ee**1.5)
     1     + ((ecm**3)/6.d0 + 5.d-1*ecm + 2.5d-1/ecm)
     2       /dsqrt(ee)
      endif
c
 900  return
      end
c
c
c
c
      double precision function I_pm(ecm, eps1, eps2, c)
      implicit none
      double precision ecm, eps1, eps2, c
c
      double precision d2
c
      if (c.gt.1.d-40) then
         d2 = eps1*eps2 + c*(ecm**2)
         I_pm = dlog(ecm*dsqrt(c) + dsqrt(d2))/dsqrt(c)
      else if (c.lt.-1.d-40) then
         d2 = -c/(eps1*eps2)
         I_pm = asin(ecm*dsqrt(d2))/dsqrt(-c)
      else
         I_pm = 0.d0
      endif
c
 900  return
      end
c
c
c      This subroutine calculates
c      the pair annihilation rates
c
c
      subroutine pa_calc(j, k, ne_local)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer j, k
      double precision ne_local
c
c
      integer i
c
      double precision gamma, pa_el, pa_pos,
     3                 Delta_ne, Delta_np, Delta_g
c
c
      if (k.eq.1) then
         if (j.eq.1) open(17, file='pa_1.dat', status='unknown')
         if (j.eq.2) open(18, file='pa_2.dat', status='unknown')
      endif
c
      Delta_ne = 0.d0
      Delta_np = 0.d0
      do 100 i = 1, num_nt
         gamma = 1.d0 + gnt(i)
         if (f_nt(j, k, i).gt.0.d0) then
            dne_pa(j, k, i) = -ne_local*f_nt(j, k, i)*pa_el(j, k, gamma)
         else
            dne_pa(j, k, i) = 0.d0
         endif
         if (n_pos(j, k, i).gt.0.d0) then
            dnp_pa(j, k, i) = -n_pos(j, k, i)*ne_local*pa_pos(j,k,gamma)
         else
            dnp_pa(j, k, i) = 0.d0
         endif
c
         if (k.gt.1) goto 50
         if (j.eq.1) then
            write(17, 200) gamma, dmax1(1.d-20, -dne_pa(j, k, i)), 
     1                     dmax1(1.d-20, -dnp_pa(j, k, i))
         else if (j.eq.2) then
            write(18, 200) gamma, dmax1(1.d-20, -dne_pa(j, k, i)), 
     1                     dmax1(1.d-20, -dnp_pa(j, k, i))
         endif
 50      continue
c
         if (i.lt.num_nt) then
            Delta_g = gnt(i+1) - gnt(i)
            Delta_ne = Delta_ne + dne_pa(j, k, i)*Delta_g
            Delta_np = Delta_np + dnp_pa(j, k, i)*Delta_g
         endif
 100  continue
c
 200  format (e12.3,1x,e12.3,1x,e12.3)
c
      if (k.eq.1) then
         if (j.eq.1) then
            close(17)
         else if (j.eq.2) then
            close(18)
         endif
      endif
c
      return
      end
c
c
c
c
      double precision function pa_el(j, k, ge)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer j, k
      double precision ge
c
c
      integer i
c
      double precision gp
      double precision Delta_g, sum, vsigma
c
c
      sum = 0.d0
      do 100 i = 1, num_nt-1
         if (n_pos(j, k, i).gt.0.d0) then
            gp = gnt(i) + 1.d0
            Delta_g = gnt(i+1) - gnt(i)
            sum = sum + n_pos(j, k, i)*Delta_g*vsigma(ge, gp)
         endif
 100  continue
      pa_el = sum
c
      return
      end
c
c
c
c
      double precision function pa_pos(j, k, gp)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer j, k
      double precision gp
c
c
      integer i
c
      double precision ge
      double precision Delta_g, sum, vsigma
c
c
      sum = 0.d0
      do 100 i = 1, num_nt-1
         if (f_nt(j, k, i).gt.0.d0) then
            ge = gnt(i) + 1.d0
            Delta_g = gnt(i+1) - gnt(i)
            sum = sum + f_nt(j, k, i)*Delta_g*vsigma(ge, gp)
         endif
 100  continue
      pa_pos = sum
c
      return
      end
c
c
c
      double precision function vsigma(ge, gp)
      implicit none
      double precision ge, gp
c
      double precision gcm_min, gcm_max, f_vs, be, bp
      double precision gmin2, gmax2
c
      be = dsqrt(1.d0 - 1.d0/(ge**2))
      bp = dsqrt(1.d0 - 1.d0/(gp**2))
      gmin2 = 5.d-1*(1.d0 + ge*gp*(1.d0 - be*bp))
      if (gmin2.gt.1.00002d0) then
         gcm_min = dsqrt(gmin2)
      else
         gcm_min = 1.00001d0
      endif
      gmax2 = 5.d-1*(1.d0 + ge*gp*(1.d0 + be*bp))
      if (gmax2.gt.1.00002d0) then
         gcm_max = dsqrt(gmax2)
      else
         gcm_max = 1.00001d0
      endif
c
      if (gcm_max.gt.gcm_min) then
         vsigma = 7.48d-15*(f_vs(gcm_max) - f_vs(gcm_min))
     1           /(be*bp*((ge*gp)**2))
      else
         vsigma = 0.d0
      endif
c
      return
      end
c
c
c
      double precision function f_vs(gcm)
      double precision gcm
c      
      double precision bcm, L
c
      bcm = dsqrt(1.d0 - 1.d0/(gcm**2))
      L = dlog((1.d0 + bcm)/(1.d0 - bcm))
c
      f_vs = (bcm**3)*(gcm**2)*L - 2.d0*(gcm**2) + 7.5d-1*(L**2)
c
      return
      end
c
c
c      This subroutine smoothes the calculated
c       internal photon spectrum by fitting a
c     spectral shape n(E) = N E^(-alpha) e^(-E/E_0)
c
c
      subroutine nph_smooth(j, k, Te)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer j, k
      double precision Te

c
      integer i, k1, l, m, n1, n2
c
      double precision f_s, N, a, E0, chi2, chi2min, y
      double precision Nmin, E0min, amin, N0, E00, a0, dchi2
c
c
c
      n1 = 2
      n2 = 10
c
 100  if ((n_ph(n2, j, k).le.1.d0).or.(n_ph(n1, j, k).le.1.d0)) 
     1     goto 900
c
      a0 = dlog(n_ph(n1, j, k)/n_ph(n2, j, k))
     1    /dlog(E_gg(n2)/E_gg(n1))
      if (((a0.lt.0.d0).or.(a0.gt.4.0d0)).and.(n2.lt.15)) then
         n2 = n2 + 1
         goto 100
      endif
      if (n2.eq.15) goto 900
      if (a0.lt.1.d-2) a0 = 1.d-2
      if (a0.gt.4.0d0) a0 = 4.0d0
      N0 = n_ph(3, j, k)
      E00 = Te
c
c      write(*,*) 'a0 = ',a0
c      write(*,*) 'N0 = ',N0
c      write(*,*) 'E00 = ',E00
c
      chi2min = 1.d50
c
      N = 5.d-1*N0
      do 450 k1 = 1, 21
         a = a0 - 5.d-1
         do 400 l = 1, 13
            E0 = 3.5d-1*E00
            do 350 m = 1, 16
               chi2 = 0.d0
               do 300 i = 1, n_gg
                  y = E_gg(i)/E0
                  if (y.lt.2.d1) then
                     f_s = N/(((E_gg(i)/E_gg(3))**a)*dexp(y))
                  else
                     f_s = 0.d0
                  endif
                  if ((f_s.gt.1.d0).and.
     1                (n_ph(i, j, k).gt.1.d0)) then
                     dchi2 = ((n_ph(i, j, k) - f_s)**2)/f_s
                     chi2 = chi2 + dchi2
                  endif
 300           continue
c
               if (chi2.le.chi2min) then
                  E0min = E0
                  amin = a
                  Nmin = N
                  chi2min = chi2
               endif
               E0 = E0*1.15
 350        continue
            a = a + 5.d-2
 400     continue        
         N = N*1.075d0
 450  continue
c
      N = Nmin
      a = amin
      E0 = E0min
c
      do 500 i = 1, n_gg
         y = E_gg(i)/E0
         if (y.lt.2.d1) then
            n_ph(i, j, k) = N/(((E_gg(i)/E_gg(3))**a)*dexp(y))
         else
            n_ph(i, j, k) = 0.d0
         endif
 500   continue
c
c      write(*,*) 'Zone ',j,' N = ',N
c      write(*,*) ' a = ',a,'; E0 = ',E0
c
 900  return
      end
c
