      program main
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c     MPI variables
      logical continued
      integer status(MPI_STATUS_SIZE)
c
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      write(*,*) 'myid=',myid
      master = 0      
      inquire(file='p000_misc.dat',exist=continued)

      if(continued)then
        if(myid.eq.master)open(unit=4, file='log.txt', access='append' )
        call read_record

      else
        if(myid.eq.master) then
c
          open(unit=4, file='log.txt')
          write(4,*) 'Number of Processors:  ', numprocs
c
          call reader
c
        endif
        call setup
      endif
c
      call xec
      if(myid.eq.master) then
c
      close(4)
c
      endif
      call MPI_FINALIZE(ierr)
      write(*,*) 'myid=',myid,' end of compton2d'
c     stop ! Xuhui 7/21/08
      end
c
c

c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Fri Jun 16 16:29:29 EDT 2006
c version: 2
c Name: J. Finke
c Prints number of processors to the log file 
c at beginning of program.     
