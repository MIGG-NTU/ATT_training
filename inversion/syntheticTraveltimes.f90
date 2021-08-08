	module velmodel
        ! regular velocity model 
        integer :: nx,ny,nz
        double precision :: dx,dy,dz
        double precision,dimension(:),allocatable :: xvel,yvel,zvel
        double precision,dimension(:,:,:),allocatable :: slow
        end module
      

	
	program syntheticTraveltime
	use velmodel
!       Part A: timetable 
	double precision,dimension(:,:,:),allocatable :: timetable
        double precision :: t0,x0,y0,z0
!	arrival times
        integer :: nrec,irec
	double precision,dimension(:),allocatable :: xrec,yrec,zrec
	integer,dimension(:),allocatable :: idrec
        integer :: ievt,i,j,k,i1,j1,k1,i2,j2,k2
	integer :: ii,jj,kk,inside

	double precision :: tmp,wt(8),wx,wy,wz,tcal

	open(33,file='../data/traveltimeReceiverGathers_synthetic')

	call slowness
	allocate(timetable(nx,ny,nz))
	open(22,file='../data/traveltimeReceiverGathers')
  	do while((.true.))
	   read(22,*,end=666) ievt,z0,x0,y0,nrec 
	   write(33, '(I9,3F16.5,I12)') ievt,z0,x0,y0,nrec    
           t0 = 0.0 
           allocate(xrec(nrec),yrec(nrec),zrec(nrec),idrec(nrec))
           do i=1,nrec
                 read(22,*) ii,jj,idrec(i),zrec(i),xrec(i),yrec(i),tmp
           end do
		   print*, 'start forward modelling ', ievt 
           call forward(x0,y0,z0,timetable)
           do irec=1,nrec
             call locate_within_cube(xrec(irec),yrec(irec),zrec(irec),&
                  nx,ny,nz,xvel,yvel,zvel,i1,j1,k1,wx,wy,wz,wt,inside)
             if(inside.lt.1)then
	    	    print*,'I cannot locate the receiver in the domain: ', &
                      xrec(irec),yrec(irec),zrec(irec)
	       		stop
	     	 end if
             i2 = i1 + 1
             j2 = j1 + 1
             k2 = k1 + 1
             !print*,i1,j1,k1,i2,j2,k2
             tcal = wt(1)*timetable(i1,j1,k1)+wt(2)*timetable(i1,j2,k1) + &
                    wt(3)*timetable(i2,j2,k1)+wt(4)*timetable(i2,j1,k1) + &
                    wt(5)*timetable(i1,j1,k2)+wt(6)*timetable(i1,j2,k2) + &
                    wt(7)*timetable(i2,j2,k2)+wt(8)*timetable(i2,j1,k2)
	         write(33,100) ievt,irec,idrec(irec),zrec(irec),xrec(irec),&
        	        yrec(irec),tcal
           end do
	   	deallocate(xrec,yrec,zrec,idrec)
        end do
666     continue
        close(22)
		close(33)
        deallocate(slow,timetable)
100     FORMAT(3I9,3F16.5,F12.5)
        stop    
        end


	subroutine slowness
	use velmodel
	integer :: i,j,k
	double precision :: tmp
	open(10,file='../model/velocity3d_true') 
        read(10,*) nx,ny,nz
        allocate(xvel(nx),yvel(ny),zvel(nz),slow(nx,ny,nz))
        do j=1,ny
        do i=1,nx
        do k=1,nz
           read(10,*) xvel(i),yvel(j),zvel(k),tmp
	   slow(i,j,k) = 1.0/tmp
        end do
        end do
        end do
        close(10)
	dx = xvel(2)-xvel(1)
        dy = yvel(2)-yvel(1)
        dz = zvel(2)-zvel(1)
	return
	end

	

	subroutine forward(x0,y0,z0,timetable)
	use velmodel
	include 'parameter_forward'
	double precision :: x0,y0,z0,gap
	double precision :: timetable(nx,ny,nz),source(nx,ny,nz)
	integer :: iteration
	timetable  = valinfty
        source 	 = 0.0
	call sourceInitiation(x0,y0,z0,timetable,source,ixyz)
	iteration = 0
	gap       = 2.0*eps
        do while((iteration.lt.iterationmax).and.(gap.gt.eps))
          call fastsweep(1,nx,1,ny,1,nz,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
          call fastsweep(nx,1,1,ny,1,nz,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
          call fastsweep(1,nx,ny,1,1,nz,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
          call fastsweep(nx,1,ny,1,1,nz,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
          call fastsweep(1,nx,1,ny,nz,1,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
          call fastsweep(nx,1,1,ny,nz,1,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
          call fastsweep(1,nx,ny,1,nz,1,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
          call fastsweep(nx,1,ny,1,nz,1,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
	  iteration = iteration + 1
          !print*,'Iteration ', iteration, ' and gap is ',gap
	end do  
	return
	end

	subroutine sourceInitiation(x0,y0,z0,timetable,source,ixyz)
	use velmodel
	integer :: ixyz
	double precision :: timetable(nx,ny,nz),source(nx,ny,nz)
	double precision :: x0,y0,z0,xtmp,ytmp,ztmp
	double precision :: wx,wy,wz,wt(8)
	integer :: iloc,ii1,ii2,jj1,jj2,i1,j1,k1,i,j,k
	call locate_within_cube(x0,y0,z0,nx,ny,nz,xvel,yvel,zvel,i1,j1,k1,wx,wy,wz,wt,iloc)
        if(iloc.lt.1)then
	       print*,'I cannot locate the source in the domain: ', x0, y0, z0
	       stop
	end if
	if(wx.gt.0.5)then
	   i1 = i1 + 1
	end if
        if(wy.gt.0.5)then
	   j1 = j1 + 1
	end if
        if(wz.gt.0.5)then
	   k1 = k1 + 1
	end if
	ii1 = max(1,i1-ixyz)
	ii2 = min(nx,i1+ixyz)
        jj1 = max(1,j1-ixyz)
	jj2 = min(ny,j1+ixyz)
        kk1 = max(1,k1-ixyz)
	kk2 = min(nz,k1+ixyz)
	do i = ii1,ii2
	   xtmp = xvel(i)
	   do j = jj1,jj2
	      ytmp = yvel(j)
	      do k = kk1,kk2
		 ztmp = zvel(k)
		 source(i,j,k) = 1.0
		 timetable(i,j,k)  = sqrt((xtmp-x0)**2+(ytmp-y0)**2+(ztmp-z0)**2)*slow(i1,j1,k1)
	      end do
	   end do
	end do
	return
	end

        subroutine fastsweep(ix1,ix2,iy1,iy2,iz1,iz2,nx,ny,nz,dx,dy,dz,slow,timetable,source,gap)
	integer :: i,j,k,l
	integer :: nx,ny,nz,ix1,ix2,ix,iy1,iy2,iy,iz1,iz2,iz
	double precision :: dx,dy,dz,timetable(nx,ny,nz),slow(nx,ny,nz),source(nx,ny,nz)
        double precision :: aa,bb,cc,ss,a,b,c,delta,tnow(7)
	double precision :: tmp,tmp1,tmp2,gap
	ix = 1
	iy = 1
        iz = 1
	if(ix1.gt.ix2)then
	  ix=-1
	end if
	if(iy1.gt.iy2)then
	  iy=-1
	end if
	if(iz1.gt.iz2)then
	  iz=-1
	end if
	gap = 0.0
        do i=ix1,ix2,ix
        do j=iy1,iy2,iy
        do k=iz1,iz2,iz
	   if(source(i,j,k).lt.0.5)then
	     ss = slow(i,j,k)*slow(i,j,k)
	     if(i.eq.1)then
	       aa = timetable(i+1,j,k)
	     else if(i.eq.nx)then
	       aa = timetable(i-1,j,k)
	     else
	       aa  = min(timetable(i-1,j,k),timetable(i+1,j,k)) 
	     end if  
    	     if(j.eq.1)then
	       bb = timetable(i,j+1,k)
	     else if(j.eq.ny)then
	       bb = timetable(i,j-1,k)
	     else
	       bb  = min(timetable(i,j-1,k),timetable(i,j+1,k)) 
	     end if  
    	     if(k.eq.1)then
	       cc = timetable(i,j,k+1)
	     else if(k.eq.nz)then
	       cc = timetable(i,j,k-1)
	     else
	       cc  = min(timetable(i,j,k-1),timetable(i,j,k+1)) 
	     end if  
	     !print*,i,j,k,aa,bb,cc
	     tnow = -1.0
	     ! case 1
	     a   = 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz)
	     b   = -2.0*aa/(dx*dx) -2.0*bb/(dy*dy) -2.0*cc/(dz*dz)
	     c   = aa*aa/(dx*dx) + bb*bb/(dy*dy) + cc*cc/(dz*dz) - ss 
             delta = b*b - 4*a*c
	     if(delta.ge.0.0)then
		tmp1 = (-1.0*b-sqrt(delta))/(2.0*a)
                tmp2 = (-1.0*b+sqrt(delta))/(2.0*a)
		if((tmp2.ge.aa).and.(tmp2.ge.bb).and.(tmp2.ge.cc))then
		    tnow(1) = tmp2
		end if
		if((tmp1.ge.aa).and.(tmp1.ge.bb).and.(tmp1.ge.cc))then
		    tnow(1) = tmp1
		end if
	     end if	
             ! case 2
	     a   = 1.0/(dx*dx) + 1.0/(dy*dy) 
	     b   = -2.0*aa/(dx*dx) -2.0*bb/(dy*dy) 
	     c   = aa*aa/(dx*dx) + bb*bb/(dy*dy) - ss 
             delta = b*b - 4*a*c
	     if(delta.ge.0.0)then
                tmp1 = (-1.0*b-sqrt(delta))/(2.0*a)
                tmp2 = (-1.0*b+sqrt(delta))/(2.0*a)
		if((tmp2.ge.aa).and.(tmp2.ge.bb).and.(tmp2.le.cc))then
		    tnow(2) = tmp2
		end if
		if((tmp1.ge.aa).and.(tmp1.ge.bb).and.(tmp1.le.cc))then
		    tnow(2) = tmp1
		end if
	     end if
	     ! case 3
	     a   = 1.0/(dx*dx) + 1.0/(dz*dz)
	     b   = -2.0*aa/(dx*dx)-2.0*cc/(dz*dz)
	     c   = aa*aa/(dx*dx) + cc*cc/(dz*dz) - ss 
             delta = b*b - 4*a*c
	     if(delta.ge.0.0)then
		tmp1 = (-1.0*b-sqrt(delta))/(2.0*a)
                tmp2 = (-1.0*b+sqrt(delta))/(2.0*a)
		if((tmp2.ge.aa).and.(tmp2.le.bb).and.(tmp2.ge.cc))then
		    tnow(3) = tmp2
		end if
		if((tmp1.ge.aa).and.(tmp1.le.bb).and.(tmp1.ge.cc))then
		    tnow(3) = tmp1
		end if
	     end if   
   	     ! case 4
	     a   = 1.0/(dy*dy) + 1.0/(dz*dz)
	     b   = -2.0*bb/(dy*dy) -2.0*cc/(dz*dz)
	     c   = bb*bb/(dy*dy) + cc*cc/(dz*dz) - ss 
             delta = b*b - 4*a*c
	     if(delta.ge.0.0)then
		tmp1 = (-1.0*b-sqrt(delta))/(2.0*a)
                tmp2 = (-1.0*b+sqrt(delta))/(2.0*a)
		if((tmp2.le.aa).and.(tmp2.ge.bb).and.(tmp2.ge.cc))then
		    tnow(4) = tmp2
		end if
		if((tmp1.le.aa).and.(tmp1.ge.bb).and.(tmp1.ge.cc))then
		    tnow(4) = tmp1
		end if
	     end if
             ! case 5
	     a   = 1.0/(dx*dx) 
	     b   = -2.0*aa/(dx*dx) 
	     c   = aa*aa/(dx*dx) - ss
             delta = b*b - 4*a*c
	     if(delta.ge.0.0)then
		tmp1 = (-1.0*b-sqrt(delta))/(2.0*a)
                tmp2 = (-1.0*b+sqrt(delta))/(2.0*a)
		if((tmp2.ge.aa).and.(tmp2.le.bb).and.(tmp2.le.cc))then
		    tnow(5) = tmp2
		end if
		if((tmp1.ge.aa).and.(tmp1.le.bb).and.(tmp1.le.cc))then
		    tnow(5) = tmp1
		end if
	     end if
   	     ! case 6
	     a   =  1.0/(dy*dy) 
	     b   = -2.0*bb/(dy*dy)
	     c   = bb*bb/(dy*dy) - ss 
             delta = b*b - 4*a*c
	     if(delta.ge.0.0)then
		tmp1 = (-1.0*b-sqrt(delta))/(2.0*a)
                tmp2 = (-1.0*b+sqrt(delta))/(2.0*a)
		if((tmp2.le.aa).and.(tmp2.ge.bb).and.(tmp2.le.cc))then
		    tnow(6) = tmp2
		end if
		if((tmp1.le.aa).and.(tmp1.ge.bb).and.(tmp1.le.cc))then
		    tnow(6) = tmp1
		end if
	     end if
     	     ! case 7
	     a   = 1.0/(dz*dz)
	     b   = -2.0*cc/(dz*dz)
	     c   = cc*cc/(dz*dz) - ss 
             delta = b*b - 4*a*c
	     if(delta.ge.0.0)then
		tmp1 = (-1.0*b-sqrt(delta))/(2.0*a)
                tmp2 = (-1.0*b+sqrt(delta))/(2.0*a)
		if((tmp2.le.aa).and.(tmp2.le.bb).and.(tmp2.ge.cc))then
		    tnow(7) = tmp2
		end if
		if((tmp1.le.aa).and.(tmp1.le.bb).and.(tmp1.ge.cc))then
		    tnow(7) = tmp1
		end if
	     end if
	     tmp = timetable(i,j,k)
	     do l = 1,7	
		if((tnow(l).ge.0.0).and.(tnow(l).le.tmp))then
		    tmp = tnow(l)
		end if
	     end do
	     gap = gap + (timetable(i,j,k)-tmp)**2
	     timetable(i,j,k) = tmp
	  else
	    !print*,'here is the source ', i,j,k,timetable(i,j,k)
	  end if
        end do
        end do
        end do
	gap = sqrt(gap*dx*dy*dz)
	!print*,gap
	return
	end
                   
	
	

        subroutine locate_within_cube(xx,yy,zz,nx,ny,nz,xvel,yvel,zvel,i1,j1,k1,wx,wy,wz,wt,iloc)
      	integer :: nx,ny,nz
	double precision :: xvel(nx),yvel(ny),zvel(nz)
        integer :: i1,i2,j1,j2,k1,k2,iloc
        double precision :: xx,yy,zz,wx,wy,wz,wt(8),val1
        logical :: notfindx,notfindy,notfindz
        iloc = 0
        i1 = 1
        j1 = 1
        k1 = 1
        notfindx =.true.
        do while(notfindx.and.i1.gt.0.and.i1.lt.nx)
          i2 = i1+1
          if(xx.lt.xvel(i1))then
          i1 = i1-1
          else
            if(xx.le.xvel(i2))then
            notfindx =.false.
            wx = (xx-xvel(i1))/(xvel(i2)-xvel(i1))
            else
            i1 = i1+1
            end if
          end if
        end do
        notfindy =.true.
        do while(notfindy.and.j1.gt.0.and.j1.lt.ny)
          j2 = j1+1
          if(yy.lt.yvel(j1))then
           j1 = j1-1
          else
            if(yy.le.yvel(j2))then
            notfindy =.false.
            wy = (yy-yvel(j1))/(yvel(j2)-yvel(j1))
            else
            j1 = j1+1
            end if
          end if
        end do
        notfindz =.true.
        do while(notfindz.and.k1.gt.0.and.k1.lt.nz)
          k2 = k1+1
          if(zz.lt.zvel(k1))then
            k1 = k1-1
          else
            if(zz.le.zvel(k2))then
            notfindz =.false.
            wz = (zz-zvel(k1))/(zvel(k2)-zvel(k1))
            else
            k1 = k1+1
            end if
          end if
        end do
        if((.not.notfindx).and.(.not.notfindy).and.(.not.notfindz))then
         iloc =1
         wt(1)=(1.0-wx)*(1.0-wy)*(1.0-wz)
         wt(2)=(1.0-wx)*wy*(1.0-wz)
         wt(3)=wx*wy*(1.0-wz)
         wt(4)=wx*(1.0-wy)*(1.0-wz)
         wt(5)=(1.0-wx)*(1.0-wy)*wz
         wt(6)=(1.0-wx)*wy*wz
         wt(7)=wx*wy*wz
         wt(8)=wx*(1.0-wy)*wz
        end if
        return
        end
