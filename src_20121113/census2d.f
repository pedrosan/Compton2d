      subroutine write_cens(entries, nunit_c)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer entries, nunit_c
c
c      integer ucens, iucens, ducens
c      parameter (ucens = 1000)
c      parameter(iucens = 6000)
c      parameter(ducens = 6000)  
c      parameter (ucens = 50000000)
c      parameter(iucens = 300000000)
c      parameter(ducens = 300000000) ! Xuhui ucens

c
      integer i, n, lwad, lwai, write_flag
c
c      
c
    
          lwad = 0
          lwai = 0   
          do 100 n = 1, entries
             write(nunit_c, 200) (dbufout(lwad+i), i=1,6)
             write(nunit_c, 210) (ibufout(lwai+i), i=1,6)
             lwad = lwad + 6
             lwai = lwai + 6
 100      continue
c
c 200  format (5e14.7)
c 210  format (5i5)
200   format (6e14.7)
210   format (6i5) !Xuhui
c
      return
      end
c
c
c
      subroutine read_cens(entries,nunit_c)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
c
      integer entries, nunit_c
c
c      integer ucens, iucens, ducens,
c      parameter (ucens = 1000)
c      parameter (iucens = 6000)
c      parameter (ducens = 6000)
c      parameter (ucens = 50000000)
c      parameter (iucens = 300000000)
c      parameter (ducens = 300000000) ! Xuhui ucens

c
      integer i, n, lwad, lwai
c
c
c
c      write(*,*) 'in read_cens'
          lwad = 0
          lwai = 0
          do 100 n = 1, entries
             read(nunit_c, 200) (dbufout(lwad+i), i=1,6)
             read(nunit_c, 210) (ibufout(lwai+i), i=1,6)
         lwad = lwad + 6
         lwai = lwai + 6
 100      continue
c
c 200  format (5e14.7)
c 210  format (5i5)
200   format (6e14.7)
210   format (6i5)  !Xuhui
c
      return
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 14:50:09 EDT 2006
c version: 2
c Name: J. Finke
c Removed variable nmax which is not used.  
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Fri Jun 16 12:24:11 EDT 2006
c version: 3
c Name: J. Finke
c fibran seed is now stored in census files. 
c imcfield is now determinable.     
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c Tue Sep 23 15:47 CDT 2008
c Name: Xuhui Chen
c When read_flag and write_flag is 0, census file
c read and write is turned off. The information
c is just stored in buffer.
