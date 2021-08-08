
       program local_xyz

       implicit none
       
       integer,parameter :: nevt=47,nsta=25
       integer :: ievt,ista,jevt,jsta 
       !earthquake information
       integer :: iyear,imon,iday,ihour,imin,idevt
       double precision :: sec,xlat,xlon,xdep,xmag
       !station information
       double precision :: ylat,ylon,yhgt,epi,atime
       character (len=1) :: phase
       character (len=3) :: cindex
       character (len=9) :: staname

      open(10,file='dataInFormat_D_xyz')       

      iyear = 2020 - nevt
      imon  = 4
      iday  = 9
      ihour = 8
      imin  = 8
      sec   = 8.8

      xlat  = 0.0
      xdep  = -15.0   ! in km. negative: below the sea level; postive: above the sea level  
      xmag  = 3.0
      
!     Seismic stations

      phase = 'P'
      ylat  = xlat  
      yhgt  = 0.0    ! in km. negative: below the sea level; postive: above the sea level 
      atime = 99.0

      do ievt = 1,nevt
	 iyear = iyear + 1 
	 xlon  = 5.0 + (ievt-1)*5.0
	 idevt = 100 + ievt
         jevt = ievt
         write(10,200) jevt,iyear,imon,iday,ihour,imin,sec,&
                                      xlon,xlat,xdep,xmag,nsta,idevt 
         do ista = 1,nsta
	    write(cindex,'(I3.3)') ista
	    staname='STA'//cindex
            ylon = 0.0 + (ista-1)*10.0
            epi  = sqrt((xlat-ylat)**2+(xlon-ylon)**2)
            jsta = ista
            write(10,400) jevt,jsta,staname,ylon,ylat,yhgt,phase,epi,atime
	 end do            
	  !end if
      end do
 
      close(10) 
     

200   format(I10,I5,I3,I3,I3,I3,2X,F6.3,F10.4,F12.4,F7.2,F6.2,I5,I10) 
400   format(2I7,1X,A10,3F12.4,2X,A1,F8.2,F8.3)
      stop 
      end program local_xyz 
