      program dataSelect
      implicit none
!     station part
      double precision :: qwest,qeast,qsouth,qnorth 
!     earthquake part
      double precision :: pwest,peast,psouth,pnorth,ptop,pbot,amag,bmag,cmag 
      integer :: nevtin
      
      call researchRegion(pwest,peast,psouth,pnorth,ptop,pbot,&
                      amag,bmag,cmag,qwest,qeast,qsouth,qnorth)
        !print*,pwest,peast,psouth,pnorth,ptop,pbot,qwest,qeast,qsouth,qnorth
      call earthquakesWithRecordsInResearchRegion(nevtin,pwest,peast,&
           psouth,pnorth,ptop,pbot,amag,bmag,cmag,qwest,qeast,qsouth,qnorth)
      if(nevtin.gt.0)then
      call filter(pwest,peast,psouth,pnorth,ptop,pbot,nevtin)

      call eikonalInput(nevtin)
      else
        print*,'no earthquaks in your research region'
      end if

      stop
      end


      subroutine eikonalInput(nevtin)
      integer,parameter :: nmax = 100000
      character (len=10) :: station(nmax) 
      double precision,dimension(nmax) :: ystalon,ystalat,ystahgt
      !earthquake information
      integer :: ievt,iyear,imon,iday,ihour,imin,idevt
      double precision :: sec,xlat,xlon,xdep,xmag,tori
      !station information
      double precision :: ylat,ylon,yhgt,epi,atime
      integer :: ista,iphase
      character (len=10) :: staname
      ! temporal parameters
      integer :: i,j,k,iglob,nevtin,nray
      integer,dimension(:,:),allocatable :: isrcrec
      double precision,dimension(:,:),allocatable :: tsrcrec
      double precision :: val,eps,xyzevt(nevtin,3)
      double precision :: xmin,xmax,ymin,ymax,zmin,zmax
      eps = 1.0e-4
      iglob = 0
      nray  = 0
      xmin = 1.0e+6
      xmax = -1.0*xmin
      ymin = 1.0e+6
      ymax = -1.0*ymin  
      zmin = 1.0e+6
      zmax = -1.0*zmin
      open(30,file='selectedData')
      open(31,file='sources')
      open(32,file='receivers')
      open(33,file='src-rec')
      do while(.true.)
           read(30,200,end=333) ievt,iyear,imon,iday,ihour,imin,sec,&
                              xlon,xlat,xdep,xmag,ista,idevt
	   tori = 0.0
           write(31,'(I9,4F14.5)') ievt,xdep,xlon,xlat,tori
           xyzevt(ievt,1) = xlon
	   xyzevt(ievt,2) = xlat
           xyzevt(ievt,3) = xdep
           if(xlon.gt.xmax) xmax = xlon
           if(xlon.lt.xmin) xmin = xlon
           if(xlat.gt.ymax) ymax = xlat
           if(xlat.lt.ymin) ymin = xlat
           if(xdep.gt.zmax) zmax = xdep
           if(xdep.lt.zmin) zmin = xdep
           do i = 1,ista
              read(30,401) k,staname,ylon,ylat,yhgt,iphase,epi,atime
              do j = 1,iglob
		if(station(j).eq.staname)then
                  val = 100.0*abs(ystalon(j)-ylon)&
                       +100.0*abs(ystalat(j)-ylat)&
                       +0.001*abs(ystahgt(j)-yhgt)  ! here it is in kilometer
                  if(val.gt.eps)then
                    !print*,'Warning: relocated and considered as a new station 1',J,ijk,staname,val,stalat,stalon,stahgt
                    !print*,'Warning: relocated and considered as a new station 2',J,ijk,staname,val,xstalat(ijk),xstalon(ijk),xstahgt(ijk)
                  else 
                    !print*,'same name, same place',j,station(j),staname
                    go to 900
                  end if
                end if          ! same name
              end do
              iglob = iglob + 1
              j = iglob
              station(j) = staname
              ystalon(j) = ylon
	      ystalat(j) = ylat
              ystahgt(j) = yhgt
900	      continue
              write(33,402) ievt,i,j,staname,ylon,ylat,yhgt,iphase,epi,atime
    	        if(ylon.gt.xmax) xmax = ylon
              if(ylon.lt.xmin) xmin = ylon
              if(ylat.gt.ymax) ymax = ylat
              if(ylat.lt.ymin) ymin = ylat
              if(yhgt.gt.zmax) zmax = yhgt
              if(yhgt.lt.zmin) zmin = yhgt
           end do
      end do
333   continue
      do i = 1,iglob
	write(32,'(I9,3F14.5)') i,ystahgt(i),ystalon(i),ystalat(i)
      end do
      close(30)
      close(31)
      close(32)
      close(33)
      open(15,file='evtstaminmax') 
      write(15,'(A35,F12.5,A7,F12.5)') &
             "double precision,parameter :: xmin=",xmin,", xmax=",xmax
      write(15,'(A35,F12.5,A7,F12.5)') &
             "double precision,parameter :: ymin=",ymin,", ymax=",ymax
      write(15,'(A35,F12.5,A7,F12.5)') &
             "double precision,parameter :: zmin=",zmin,", zmax=",zmax
      close(15)
      allocate(isrcrec(ievt,iglob),tsrcrec(ievt,iglob))
      isrcrec = 0
      tsrcrec = 0.0
      open(33,file='src-rec')
      do while(.true.)
	read(33,402,end=222) k,i,j,staname,ylon,ylat,yhgt,iphase,epi,atime
        if(isrcrec(k,j).gt.0)then
           isrcrec(k,j) =  isrcrec(k,j) + 1
           tsrcrec(k,j) =  tsrcrec(k,j) + atime
           print*,'more than one record', k,j
        else
           isrcrec(k,j) = 1          
           tsrcrec(k,j) = atime
        end if
      end do
222   continue
      close(33)
      open(35,file='obstime-src-gather')
      do i = 1,ievt
	 l = 0
      do j = 1,iglob
         if(isrcrec(i,j).gt.0)then
           l = l + 1
	   write(35,*) i,l,j,ystahgt(j),ystalon(j),&
                 ystalat(j),tsrcrec(i,j)/isrcrec(i,j)
         end if
      end do
      end do
      close(35)

      nray = 0
      open(35,file='obstime-rec-gather')
      do i = 1,iglob
         l = 0
      do j = 1,ievt
         if(isrcrec(j,i).gt.0)then
           l = l + 1
           nray = nray + 1
           write(35,*) i,l,j,xyzevt(j,3),xyzevt(j,1),&
                       xyzevt(j,2),tsrcrec(j,i)/isrcrec(j,i)
         end if
      end do
      end do
      close(35)

      open(43,file='nevtstaray')
      write(43,'(A26,I12)') "integer,parameter :: nevt=",ievt
      write(43,'(A26,I12)') "integer,parameter :: nsta=",iglob
      write(43,'(A26,I12)') "integer,parameter :: nray=",nray
      close(43)
      print*,'nsrc,nrec,ndata',ievt,iglob,nray
401   format(I7,A10,3F12.4,2X,I1,F8.2,F8.3)
402   format(3I7,A10,3F12.4,2X,I1,F8.2,F8.3)
200   format(I10,I5,I3,I3,I3,I3,2X,F6.3,F10.4,F12.4,F7.2,F6.2,I5,I10)
      return
      end

      subroutine filter(pwest,peast,psouth,pnorth,ptop,pbot,nevtin)
      !research region information
      integer :: nevtin 
      double precision :: pwest,peast,psouth,pnorth,ptop,pbot
      !earthquake information
      integer :: ievt,iyear,imon,iday,ihour,imin,idevt
      double precision :: sec,xlat,xlon,xdep,xmag
      !station information
      double precision :: ylat,ylon,yhgt,epi,atime
      integer :: ns,iphase
      character (len=10) :: staname
      !other parameters
      integer :: i,j,k
      !boxes
      integer :: nlat,nlon,ndep,ilat,ilon,idep
      integer,dimension(:,:,:),allocatable :: idquake,maxrec
      integer :: iselected(nevtin)
      iselected = 0
      open(10,file='../0_commandCenter/blockSelection')
      read(10,'(3F12.5)') dlat,dlon,ddep
      read(10,'(I6)') nbd
      close(10)
      nlat = max(ifix(real((pnorth - psouth)/dlat))+1,2)
      nlon = max(ifix(real((peast  - pwest)/dlon))+1,2)
      ndep = max(ifix(real((ptop   - pbot)/ddep))+1,2)
      allocate(idquake(nlat,nlon,ndep),maxrec(nlat,nlon,ndep))
      idquake = 0
      maxrec  = 0
      iselected = 0
      open(21,file='earthquakesOnlyInResearchRegion')
      do while(.true.)
         read(21,200,end=666) ievt,iyear,imon,iday,ihour,imin,sec,&
                                     xlon,xlat,xdep,xmag,ista,idevt
	 ilat = ifix(real((xlat - psouth)/dlat))+1
      	 ilon = ifix(real((xlon - pwest)/dlon))+1
         idep = ifix(real((xdep - pbot)/ddep))+1
	 if(maxrec(ilat,ilon,idep).lt.ista)then   ! find one earthquake with more recordings in the box
           j = idquake(ilat,ilon,idep)            ! get the id of the previous earthquake that has the maximum number of recordings
           if(j.gt.0)then
	     iselected(j) = 0                     ! now we do not want to select it anymore
             print*,'Earthquake ',j,' is replaced by ',ievt
	   end if
	   maxrec(ilat,ilon,idep)  = ista         
	   idquake(ilat,ilon,idep) = ievt
           iselected(ievt) = 1
         else
           print*,'a bigger one came earlier',idquake(ilat,ilon,idep)
           write(6,200) ievt,iyear,imon,iday,ihour,imin,sec,&
                               xlon,xlat,xdep,xmag,ista,idevt
	 end if
      end do
666   continue
      close(21)
      open(20,file='phaseDataInResearchRegion')
      open(30,file='selectedData')
      k = 0
      do while(.true.)
         read(20,200,end=999) ievt,iyear,imon,iday,ihour,imin,sec,&
                                       xlon,xlat,xdep,xmag,ista,idevt
	 if((iselected(ievt).gt.0).and.(ista.ge.nbd))then
           k = k + 1
           write(30,200) k,iyear,imon,iday,ihour,imin,sec,&
                              xlon,xlat,xdep,xmag,ista,idevt
           do i = 1,ista
             read(20,401) j,staname,ylon,ylat,yhgt,iphase,epi,atime
             write(30,401) j,staname,ylon,ylat,yhgt,iphase,epi,atime
           end do
	 else
           do i = 1,ista
             read(20,401) j,staname,ylon,ylat,yhgt,iphase,epi,atime
           end do
         end if
      end do
999   continue
      close(20)
      close(30)
401   format(I7,A10,3F12.4,2X,I1,F8.2,F8.3)
200   format(I10,I5,I3,I3,I3,I3,2X,F6.3,F10.4,F12.4,F7.2,F6.2,I5,I10) 
      return
      end


      subroutine earthquakesWithRecordsInResearchRegion(ievtin,pwest,peast,&
          psouth,pnorth,ptop,pbot,amag,bmag,cmag,qwest,qeast,qsouth,qnorth)
      !research region information
      include '../0_commandCenter/phasewanted'
      integer,parameter :: nsta = 100000
      double precision :: pwest,peast,psouth,pnorth,ptop,pbot,&
                       amag,bmag,cmag,qwest,qeast,qsouth,qnorth
      !earthquake information
      integer :: ievt,iyear,imon,iday,ihour,imin,idevt
      double precision :: sec,xlat,xlon,xdep,xmag
      !station information
      double precision :: ylat,ylon,yhgt,epi,atime
      character (len=1) :: phase
      integer :: ns,iphase
      character (len=10) :: staname
      !other parameters
      integer :: i,j,k,ista,ievtin,idbeg,idend
      logical :: inside
      !global station information
      character (len=10) :: station(nsta) 
      double precision,dimension(nsta) ::ystalon,ystalat,ystahgt,yepidis,yartime 
      open(10,file="../../dataInFormat_D_xyz")
      open(20,file='phaseDataInResearchRegion')
      open(21,file='earthquakesOnlyInResearchRegion')
      open(30,file='largeEarthquakes')
      open(40,file='../0_commandCenter/timePeriod')
      read(40,*) idbeg,idend
      close(40)
      ievtin = 0
      do while(.true.)
        read(10,200,end=777) ievt,iyear,imon,iday,ihour,imin,sec,&
                                      xlon,xlat,xdep,xmag,ns,idevt
        if(xmag.ge.cmag)then
           write(30,205) iyear,imon,iday,ihour,imin,sec,xlon,xlat,xdep,xmag,ns
         end if
        inside = .false.
        j = iyear*10000+imon*100+iday
	      if((j.ge.idbeg).and.(j.le.idend))then
      	if((xlon.ge.pwest).and.(xlon.le.peast))then
      	if((xlat.ge.psouth).and.(xlat.le.pnorth))then
        if((xdep.ge.pbot).and.(xdep.le.ptop))then 
        if((xmag.ge.amag).and.(xmag.le.bmag).and.(ns.ge.1))then      
	          inside = .true. 
        end if
        end if
        end if
        end if 
        end if
        if(inside)then
	  ista = 0
	  do i = 1,ns
             read(10,400) j,k,staname,ylon,ylat,yhgt,phase,epi,atime
      	     iphase = 0
      	     if(phase.eq.'P')then
        	iphase = 1
      	     else if(phase.eq.'S')then
        	iphase = 2
      	     end if
	     if(iphase.eq.ips)then
   	     if((ylon.ge.qwest).and.(ylon.le.qeast))then 
	     if((ylat.ge.qsouth).and.(ylat.le.qnorth))then    
               	ista = ista + 1
	        station(ista) = staname
      		ystalon(ista) = ylon
      		ystalat(ista) = ylat
      		ystahgt(ista) = yhgt
		yepidis(ista) = epi
		yartime(ista) = atime 	
             end if
             end if
	     end if
          end do
	  if(ista.gt.0)then
   	    ievtin = ievtin + 1
	    write(20,200) ievtin,iyear,imon,iday,ihour,imin,sec,&
                               xlon,xlat,xdep,xmag,ista,idevt
	    write(21,200) ievtin,iyear,imon,iday,ihour,imin,sec,&
                               xlon,xlat,xdep,xmag,ista,idevt
	    do i = 1,ista
	       write(20,401) i,station(i),ystalon(i),ystalat(i),&
                       ystahgt(i),ips,yepidis(i),yartime(i)
  	    end do
	  end if
	else
          do i = 1,ns
             read(10,400) j,k,staname,ylon,ylat,yhgt,phase,epi,atime
          end do
        end if
      end do
777   continue
      close(10)
      close(20)
      close(21)
      close(30)
400   format(2I7,1X,A10,3F12.4,2X,A1,F8.2,F8.3)
401   format(I7,A10,3F12.4,2X,I1,F8.2,F8.3)
200   format(I10,I5,I3,I3,I3,I3,2X,F6.3,F10.4,F12.4,F7.2,F6.2,I5,I10) 
205   format(I5,4I3,F9.3,3F12.4,F4.1,I12)
      return
      end

 

      subroutine researchRegion(pwest,peast,psouth,pnorth,ptop,pbot,&
                          amag,bmag,cmag,qwest,qeast,qsouth,qnorth)
      double precision :: pwest,peast,psouth,pnorth,ptop,pbot,&
                          amag,bmag,cmag,qwest,qeast,qsouth,qnorth
      double precision,parameter :: eps = 0.00001
      open(11,file='../0_commandCenter/earthquakeRegion')
      read(11,'(9F12.4)') pwest,peast,psouth,pnorth,pbot,ptop,&
                          amag,bmag,cmag
      close(11)
      open(12,file='../0_commandCenter/stationRegion')
      read(12,*) qwest,qeast,qsouth,qnorth
      close(12)
      pwest  = pwest  + eps
      peast  = peast  - eps
      psouth = psouth + eps
      pnorth = pnorth - eps
      ptop   = ptop + eps
      pbot   = pbot - eps
      qwest  = qwest  + eps
      qeast  = qeast  - eps
      qsouth = qsouth + eps
      qnorth = qnorth - eps
      return
      end 

