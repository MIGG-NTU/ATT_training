	module velmodel
        ! regular velocity model 
        integer :: nx,ny,nz
        double precision :: dx,dy,dz
        double precision,dimension(:),allocatable :: xvel,yvel,zvel
        double precision,dimension(:,:,:),allocatable :: slow,sref
        end module
      
	module mesh      
      	integer :: nset,ninv
      	! regular parameterization
      	integer,dimension(:),allocatable :: invx,invy,invz,istartat
      	double precision,dimension(:,:),allocatable :: xinv,yinv,zinv
        end module

	
	program adjointTraveltimeTomography
	use velmodel
	use mesh
	include "./parameter_inversion"
!       Part A: timetable and adjoint field
	double precision,dimension(:,:,:),allocatable :: timetable,adjntable,kernel
	double precision,parameter :: eps=1.0e-8
    double precision :: t0,x0,y0,z0
!	arrival times
    integer :: nrec,irec,idrec
	double precision,dimension(:),allocatable :: xrec,yrec,zrec,tobs
        integer :: ievt,i,j,k
	integer :: ii,jj,kk
	integer :: iter
	logical :: notoptimal
!       regularization 
	integer :: lambda,leta
	double precision :: valmax,dslowmax

        ! nonlinear congjugate-gradient method 
        ! gk
	double precision :: chi,step
	double precision,dimension(:),allocatable :: gk,pk

	double precision :: f1,g1,f2,gknorm
	double precision :: tmpa,tmpb,tmp
	integer :: mode
	character (len=3) :: citer
	character (len=13) :: filename,fileb
    ! This algoritm considers tobs as the time on the clock. The earthquake occurred at the time around 0. 
    ! Thus, tobs = T^o + tau^o. Note that we have no idea about T^o.

    dslowmax = updatemax

	open(66,file='objectiveFunction')

	call slowness
	allocate(timetable(nx,ny,nz),adjntable(nx,ny,nz),kernel(nx,ny,nz))

	call inversionGrid
	allocate(gk(ninv),pk(ninv))
	iter = 0
    notoptimal = .true.
	do while((iter.lt.niter).and.(notoptimal))
        if(iter.gt.0)then
            f2 = f1
	   	else
	    	f2 = 0.0
        end if
        gk  = 0.0
        chi = 0.0
        if(iter.eq.0)then
           kernel = 0.0
	   	end if
        write(citer,'(I3.3)') iter
        filename = 'residual'//citer
        open(2,file=filename)
        open(22,file='../data/traveltimeReceiverGathers')
  	    do while((.true.))
	      read(22,*,end=666) ievt,z0,x0,y0,nrec    
	      t0 = 0.0 
	      allocate(xrec(nrec),yrec(nrec),zrec(nrec),tobs(nrec))
	      do i=1,nrec
	         read(22,*) ii,jj,idrec,zrec(i),xrec(i),yrec(i),tobs(i)
	      end do 
          call forward(x0,y0,z0,timetable)
	      !print*, 'forward modelling is done'
	      mode = 2
	      call adjoint(t0,timetable,adjntable,gk,nrec,xrec,yrec,zrec,tobs,tbd,tabs,chi,mode)
	      if(iter.eq.0)then
   	        kernel = kernel + adjntable
            write(citer,'(I3.3)') ievt
	        fileb = 'adjtfield'//citer
	        open(69,file=fileb)
	        do jj=1,nz
	        do ii=1,nx
		       write(69,*) xvel(ii),yvel(6),zvel(jj),adjntable(ii,6,jj)*slow(ii,6,jj)*slow(ii,6,jj)
	      	end do
	      	end do
	      	close(69)
	      end if
	      print*, 'Iteration, Earthquake (virtual source), accumulated chi ', iter, ievt, chi
		  print*, '  '
	      deallocate(xrec,yrec,zrec,tobs)
        end do
666	   	continue
	   	close(22)
        close(2)
	   	gk = gk/nset
        pk = -1.0*gk
	    f1 = chi
	    g1 = 0.0
	    do i = 1,ninv
	       g1 = g1 + gk(i)*pk(i)
        end do
        step = -0.50*f1/g1
        if((iter.gt.0).and.(f1.gt.f2))then
            dslowmax = dslowmax*0.9
            print*,'previous and current chi', f2, f1
            print*,'maximum slowness is adjusted'
        end if
        call slownessupdate(step,pk,dslowmax,valmax)
	   	print*,'new slowness model is found'
        print*,'real & upper bound max perturbation',valmax,dslowmax
        gknorm = 0.0	
	    do i=1,ninv
	       gknorm = gknorm + gk(i)*gk(i)
	    end do
        if(gknorm.le.eps)then
	        notoptimal = .false.
        end if
	    write(66,*) iter,f1,gknorm
	    print*,'current iteration index (0, 1, ...), chi and gknorm ',iter,f1,gknorm
  	    if(iter.eq.0)then
           	open(70,file='kernel')
	       	do jj=1,nz
           	do ii=1,nx
	   		  write(70,*) xvel(ii),yvel(6),zvel(jj),kernel(ii,6,jj)*slow(ii,6,jj)*slow(ii,6,jj)
	   		end do
	   		end do
	   		close(70)
	    end if
        iter = iter + 1
	    write(citer,'(I3.3)') iter
	    filename = 'velocity3d'//citer
	    open(98,file=filename)
	    write(98,*) nx,ny,nz
	    do j=1,ny
	    do i=1,nx
 	    do k=1,nz
	      write(98,'(F16.5)') 1.0/slow(i,j,k) 
	    end do
	    end do
	    end do
	    close(98)
      	end do
		close(66)
	deallocate(slow,sref,timetable,adjntable,kernel)
	deallocate(gk,pk,xvel,yvel,zvel)
	deallocate(invx,invy,invz,istartat,xinv,yinv,zinv)
    stop    
    end


	subroutine slowness
	use velmodel
	integer :: i,j,k
	double precision :: tmp
	open(10,file='../model/velocity3d') 
        read(10,*) nx,ny,nz
        allocate(xvel(nx),yvel(ny),zvel(nz),slow(nx,ny,nz),sref(nx,ny,nz))
        do j=1,ny
        do i=1,nx
        do k=1,nz
           read(10,*) xvel(i),yvel(j),zvel(k),tmp
	   slow(i,j,k) = 1.0/tmp
        end do
        end do
        end do
        close(10)
	sref = slow
	dx = xvel(2)-xvel(1)
        dy = yvel(2)-yvel(1)
        dz = zvel(2)-zvel(1)
	return
	end

	subroutine inversionGrid
	use mesh
	integer :: iset,i,j,k,imax,jmax,kmax
	open(60,file='../mesh/multiple-grid')	
	read(60,*) nset,imax,jmax,kmax
	allocate(invx(nset),invy(nset),invz(nset),istartat(nset))
	allocate(xinv(imax,nset),yinv(jmax,nset),zinv(kmax,nset))
	ninv = 0
	do iset = 1,nset	
	   read(60,*) i,invx(iset),invy(iset),invz(iset)
	   istartat(iset) = ninv 
	   ninv = ninv+invx(iset)*invy(iset)*invz(iset)
	   do k=1,invz(iset)
           do j=1,invy(iset)
           do i=1,invx(iset)
	      read(60,*) xinv(i,iset),yinv(j,iset),zinv(k,iset)
	   end do
	   end do
	   end do
	end do
	close(60)
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
	i1 = 1
        j1 = 1
        k1 = 1
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
                   
	subroutine adjoint(torigin,timetable,adjntable,skernel,nrec,xrec,yrec,zrec,tobs,tbd,tabs,chi,mode)
	use velmodel
	use mesh
	include 'parameter_adjoint' 
        integer :: mode,nrec
	double precision :: timetable(nx,ny,nz),adjntable(nx,ny,nz),source(nx,ny,nz)
	double precision :: xrec(nrec),yrec(nrec),zrec(nrec),tobs(nrec)
	double precision :: skernel(ninv)
	integer :: i,j,k,m,n,iset,nx1,ny1,nz1
	integer :: i1,j1,k1,i2,j2,k2,inside,iteration,irec,iuse
	double precision :: wx,wy,wz,wt(8),tcal,chi
	double precision :: tmp,gap,torigin,errnum,tbd,tabs

	double precision,dimension(:),allocatable :: xinvtmp,yinvtmp,zinvtmp
	integer :: invxtmp,invytmp,invztmp
        source = 0.0   
        do irec=1,nrec
           i1 = 1
           j1 = 1
           k1 = 1
           call locate_within_cube(xrec(irec),yrec(irec),zrec(irec),&
                 nx,ny,nz,xvel,yvel,zvel,i1,j1,k1,wx,wy,wz,wt,inside)
           if(inside.lt.1)then
	       print*,'I cannot locate the receiver in the domain: ', xr, yr, zr
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
           tmp = tcal + torigin - tobs(irec)
           if((abs(tmp).le.tabs).and.(abs(tmp/tcal).le.tbd))then
	     source(i1,j1,k1) = source(i1,j1,k1) + wt(1)*tmp/(dx*dy*dz)
             source(i1,j2,k1) = source(i1,j2,k1) + wt(2)*tmp/(dx*dy*dz)
             source(i2,j2,k1) = source(i2,j2,k1) + wt(3)*tmp/(dx*dy*dz)
             source(i2,j1,k1) = source(i2,j1,k1) + wt(4)*tmp/(dx*dy*dz)
             source(i1,j1,k2) = source(i1,j1,k2) + wt(5)*tmp/(dx*dy*dz)
             source(i1,j2,k2) = source(i1,j2,k2) + wt(6)*tmp/(dx*dy*dz)
             source(i2,j2,k2) = source(i2,j2,k2) + wt(7)*tmp/(dx*dy*dz)
             source(i2,j1,k2) = source(i2,j1,k2) + wt(8)*tmp/(dx*dy*dz)
	     chi = chi + 0.5*tmp*tmp
             iuse = 1
             !print*,tcal,tobs(irec),tobs(irec)-tcal
           else
             iuse = 0
           end if
           if(mode.eq.2)then
             write(mode,*) tcal + torigin, tobs(irec),tmp,iuse
           end if
           !write(6,*) mode,i1,j1,k1,xrec(irec),yrec(irec),zrec(irec),tcal,tobs(irec)-tcal
        end do

	if(mode.eq.1)then
	  return
	end if

        adjntable = valinfty
        do j=1,ny
        do i=1,nx
        do k=1,nz
           if((j.eq.1).or.(j.eq.ny).or.(i.eq.1).or.(i.eq.nx).or.(k.eq.1).or.(k.eq.nz))then
	       adjntable(i,j,k)=0
	   end if
        end do
        end do
        end do
        nx1 = nx-1
        ny1 = ny-1
	nz1 = nz-1
  	iteration = 0
	gap = 2.0*eps
	do while((iteration.lt.iterationmax).and.(gap.gt.eps))
             call sweep(2,nx1,2,ny1,2,nz1,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
             call sweep(nx1,2,2,ny1,2,nz1,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
             call sweep(2,nx1,ny1,2,2,nz1,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
             call sweep(nx1,2,ny1,2,2,nz1,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
             call sweep(2,nx1,2,ny1,nz1,2,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
             call sweep(nx1,2,2,ny1,nz1,2,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
             call sweep(2,nx1,ny1,2,nz1,2,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
             call sweep(nx1,2,ny1,2,nz1,2,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
	     iteration = iteration + 1
             !print*,'Adjoint Iteration ', iteration, ' and gap is ',gap
	end do

	do iset=1,nset
	invxtmp = invx(iset)
	invytmp = invy(iset)
	invztmp = invz(iset)
	allocate(xinvtmp(invxtmp),yinvtmp(invytmp),zinvtmp(invztmp))
	xinvtmp = xinv(1:invxtmp,iset)
        yinvtmp = yinv(1:invytmp,iset)
	zinvtmp = zinv(1:invztmp,iset)
	i1 = 1
	j1 = 1
	k1 = 1
        do j=2,ny1
        do i=2,nx1
        do k=2,nz1
	   call locate_within_cube(xvel(i),yvel(j),zvel(k),invxtmp,invytmp,invztmp,&
                            xinvtmp,yinvtmp,zinvtmp,i1,j1,k1,wx,wy,wz,wt,inside)
	   if(inside.gt.0)then
             do n=1,8
               if(n.eq.1)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*(j1-1)+i1
               elseif(n.eq.2)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*j1+i1
               elseif(n.eq.3)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*j1+i1+1
               elseif(n.eq.4)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*(j1-1)+i1+1
               elseif(n.eq.5)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*(j1-1)+i1
               elseif(n.eq.6)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*j1+i1
               elseif(n.eq.7)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*j1+i1+1
               elseif(n.eq.8)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*(j1-1)+i1+1
               else
               end if
	       m = istartat(iset) + m
               skernel(m) = skernel(m) + adjntable(i,j,k)*slow(i,j,k)*slow(i,j,k)*dx*dy*dz*wt(n)
             end do
           else
             print*,'outside of the inversion domain ', xvel(i),yvel(j),zvel(k),i1,j1,k1
             stop
           end if        
	end do        
	end do
        end do
	deallocate(xinvtmp,yinvtmp,zinvtmp)
	end do

        return
	end 



	subroutine sweep(ix1,ix2,iy1,iy2,iz1,iz2,nx,ny,nz,dx,dy,dz,timetable,adjntable,source,gap)
        integer :: nx,ny,nz,ix1,ix2,ix,iy1,iy2,iy,iz1,iz2,iz
        double precision :: dx,dy,dz,timetable(nx,ny,nz),adjntable(nx,ny,nz),source(nx,ny,nz)
        double precision :: a1,ap1,am1,a2,ap2,am2
        double precision :: b1,bp1,bm1,b2,bp2,bm2
        double precision :: c1,cp1,cm1,c2,cp2,cm2
        double precision :: d,e,f,g,eps,tmp,gap,gapmax
        eps = 1.0e-6
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
          gapmax = 0.0
          do i=ix1,ix2,ix
          do j=iy1,iy2,iy
          do k=iz1,iz2,iz
             a1  = -1.0*(timetable(i,j,k)-timetable(i-1,j,k))/dx    ! 1 --> -1/2   ! why -T? Make sure lambda is along the negative gradient direction of the traveltime field
             ap1 = (a1 + abs(a1))/2.0
             am1 = (a1 - abs(a1))/2.0  
             a2  = -1.0*(timetable(i+1,j,k)-timetable(i,j,k))/dx    ! 2 --> +1/2
             ap2 = (a2 + abs(a2))/2.0
             am2 = (a2 - abs(a2))/2.0
             b1  = -1.0*(timetable(i,j,k)-timetable(i,j-1,k))/dy    ! 1 --> -1/2
             bp1 = (b1 + abs(b1))/2.0
             bm1 = (b1 - abs(b1))/2.0  
             b2  = -1.0*(timetable(i,j+1,k)-timetable(i,j,k))/dy    ! 2 --> +1/2
             bp2 = (b2 + abs(b2))/2.0
             bm2 = (b2 - abs(b2))/2.0
             c1  = -1.0*(timetable(i,j,k)-timetable(i,j,k-1))/dz    ! 1 --> -1/2
             cp1 = (c1 + abs(c1))/2.0
             cm1 = (c1 - abs(c1))/2.0  
             c2  = -1.0*(timetable(i,j,k+1)-timetable(i,j,k))/dz    ! 2 --> +1/2
             cp2 = (c2 + abs(c2))/2.0
             cm2 = (c2 - abs(c2))/2.0
             d   = (ap2 - am1)/dx + (bp2 - bm1)/dy + (cp2 - cm1)/dz
             if(abs(d).lt.eps)then
                adjntable(i,j,k) = 0.0
                !print*,'d value is very small ',d,i,j,k,a2,b2,c2
             else
  				e   = (ap1*adjntable(i-1,j,k)-am2*adjntable(i+1,j,k))/dx + &
                      (bp1*adjntable(i,j-1,k)-bm2*adjntable(i,j+1,k))/dy + &
                      (cp1*adjntable(i,j,k-1)-cm2*adjntable(i,j,k+1))/dz
                f   = (e + source(i,j,k))/d
                g   = adjntable(i,j,k)
                !print*,d,e,f,g
                if(g.gt.f)then
                  gap = g - f
                  adjntable(i,j,k) = f
                  if(gap.gt.gapmax)then
                     gapmax = gap
                  end if
                end if
             end if
          end do
          end do
          end do
          gap = gapmax
          !print*,'gap',gap
        return
        end

	

        subroutine locate_within_cube(xx,yy,zz,nx,ny,nz,xvel,yvel,zvel,i1,j1,k1,wx,wy,wz,wt,iloc)
      	integer :: nx,ny,nz
		double precision :: xvel(nx),yvel(ny),zvel(nz)
        integer :: i1,i2,j1,j2,k1,k2,iloc
        double precision :: xx,yy,zz,wx,wy,wz,wt(8),val1
        logical :: notfindx,notfindy,notfindz
        iloc = 0
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


	subroutine slownessupdate(step,pk,dslowmax,vmax)
	use velmodel
	use mesh
	integer :: i,j,k,iset,n,i1,j1,k1,inside
	double precision :: slowtmp(nx,ny,nz),val,vmax
	double precision :: step,pk(ninv)
	double precision :: wx,wy,wz,wt(8),dslowmax
	double precision,dimension(:),allocatable :: xinvtmp,yinvtmp,zinvtmp
	integer :: invxtmp,invytmp,invztmp
	slowtmp = 0.0
	do iset=1,nset
	i1 = 1
	j1 = 1
	k1 = 1
	invxtmp = invx(iset)
	invytmp = invy(iset)
	invztmp = invz(iset)
	allocate(xinvtmp(invxtmp),yinvtmp(invytmp),zinvtmp(invztmp))
	xinvtmp = xinv(1:invxtmp,iset)
        yinvtmp = yinv(1:invytmp,iset)
	zinvtmp = zinv(1:invztmp,iset)
        do j=1,ny
        do i=1,nx
        do k=1,nz
	   call locate_within_cube(xvel(i),yvel(j),zvel(k),invxtmp,invytmp,invztmp,&
                            xinvtmp,yinvtmp,zinvtmp,i1,j1,k1,wx,wy,wz,wt,inside)
	   if(inside.gt.0)then
	     val = 0.0
             do n=1,8
               if(n.eq.1)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*(j1-1)+i1
               elseif(n.eq.2)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*j1+i1
               elseif(n.eq.3)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*j1+i1+1
               elseif(n.eq.4)then
                   m = invx(iset)*invy(iset)*(k1-1)+invx(iset)*(j1-1)+i1+1
               elseif(n.eq.5)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*(j1-1)+i1
               elseif(n.eq.6)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*j1+i1
               elseif(n.eq.7)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*j1+i1+1
               elseif(n.eq.8)then
                   m = invx(iset)*invy(iset)*k1+invx(iset)*(j1-1)+i1+1
               else
               end if
	       m = istartat(iset) + m
               val = val + wt(n)*step*pk(m)
             end do
           else
             print*,'outside the domain: check the inversion grid'
             stop
           end if  
	   slowtmp(i,j,k) = slowtmp(i,j,k) + val      
	end do        
	end do
        end do
	deallocate(xinvtmp,yinvtmp,zinvtmp)
	end do
        slowtmp = slowtmp/nset
        vmax = -1.0
        do j=1,ny
        do i=1,nx
        do k=1,nz
           val = abs(slowtmp(i,j,k))
           if(val.gt.vmax)then
                vmax = val
           end if 
        end do
        end do
        end do
        if(vmax.gt.dslowmax)then
           step = dslowmax/vmax
        else
           step = 1.0
        end if
        do j=1,ny
        do i=1,nx
        do k=1,nz
           slowtmp(i,j,k) = 1.0 + step*slowtmp(i,j,k)
	   slow(i,j,k)  = slow(i,j,k)*slowtmp(i,j,k)
        end do
        end do
        end do
	return
	end

