    module mesh
    ! regular parameterization
    double precision,dimension(:),allocatable :: outx,outy,outz,vstartz
    double precision,dimension(:,:,:),allocatable :: vstart3d,vinvert3d,outdv3d,outray3d
    end module
      
    program vel3d
    use mesh
    include '../commandCenter/region'
    double precision,parameter :: eps = 1.0E-8
    double precision :: dxy,xbegline,ybegline,xendline,yendline,dxline,dyline
    double precision :: xx,yy,zz,dis ! coordinate transformation      
    
    integer :: i,j,k,iL,nxy
    double precision :: val1,val2,val3,val4
    double precision :: vmax,vmin,vave,dvbd 
    
    integer :: imodel
    character (len=3) :: fmtmodel,fmtline
    character (len=80) :: vfile0,vfile1
    character (len=19) :: cutfile
    
    print*,'please input model index starting from 0'
    read(*,*) imodel
    write(fmtmodel,'(I3.3)') imodel
    vfile0 = "../model/velocity3d"
    vfile1 = "../inversion/velocity3d"//fmtmodel
 
    print*,'Relative velocity perturbation bound 0.01? 0.05?'
    read(*,*) dvbd


! Read the output
   open(12,file=vfile0)
   allocate(outx(nx),outy(ny),outz(nz),vstartz(nz),vinvert3d(nx,ny,nz))
   allocate(vstart3d(nx,ny,nz),outdv3d(nx,ny,nz),outray3d(nx,ny,nz))
   vstartz = 0.0
   read(12,*) i,j,k
   do j=1,ny
      do i=1,nx
       	 do k=1,nz
       	    read(12,*) outx(i),outy(j),outz(k),vstart3d(i,j,k)
            vstartz(k) = vstartz(k) + vstart3d(i,j,k)
          end do
      end do
   end do
   close(12)
   vstartz = 1.0*vstartz/(nx*ny) 
   open(13,file='vz1d')
   do k=1,nz
      !if(vstartz(k).gt.7.8) vstartz(k) = 7.8
      write(13,*) outz(k),vstartz(k),0.95*vstartz(k),1.05*vstartz(k)
   end do
   close(13)
   open(12,file=vfile1)
   read(12,*) i,j,k
   do j=1,ny
      do i=1,nx
       	 do k=1,nz      		 
       	    read(12,*) vinvert3d(i,j,k)
            outdv3d(i,j,k) = 100.0*(vinvert3d(i,j,k)-vstart3d(i,j,k))/vstart3d(i,j,k)
            !outdv3d(i,j,k) = 100.0*(val2 - vstartz(k))/vstartz(k)
          end do
      end do
   end do
   close(12)
!************************************************ 
    open(101,file='lineEnds')    
    iL = 0
    do while(.true.)
       read(101,*,end=999) xbegline,ybegline,xendline,yendline
       iL = iL + 1
       write(fmtline,'(I3.3)') iL
       cutfile = "dvelCut"//fmtline
       open(10,file=cutfile)
 
       call epidis(xbegline,ybegline,xendline,yendline,dis)
       print*,'Line Length', dis
       dxy = dz
       nxy = ifix(real(abs(dis+0.1*dxy)/dxy))+1 
       dxline  = (xendline - xbegline)/nxy
       dyline  = (yendline - ybegline)/nxy
       vave = 0.0
       vmax = -1000.0
       vmin = 100000.0  
       i1 = 1
       j1 = 1
       k1 = 1
       do i = 1,nxy
          xx = xbegline + (i-1)*dxline
          yy = ybegline + (i-1)*dyline
          call epidis(xx,yy,xbegline,ybegline,dis)
          do k = nz,1,-1
             zz = outz(k)
             call locate_within_cube(nx,ny,nz,xx,yy,zz,i1,j1,k1,val1,val2,val3,val4,iloc)
             vave = vave + val2
             if(val2.lt.vmin) vmin = val2
             if(val2.gt.vmax) vmax = val2
             if(val2.gt.dvbd) val2 = dvbd
             if(val2.lt.(-1.0*dvbd)) val2 = -1.0*dvbd
             write(10,*) dis,-1.0*zz,val1,val2,val3,val4
          end do
       end do
       close(10)
       print*,vmin,vave/(nxy*nz),vmax       
    end do
    close(101)

999 continue
 
    deallocate(outx,outy,outz,vstartz,vstart3d,vinvert3d,outdv3d,outray3d)   
   
    stop

    end

      subroutine locate_within_cube(nx,ny,nz,xx,yy,zz,i1,j1,k1,val1,val2,val3,val4,iloc)
      use mesh
      integer :: i1,i2,j1,j2,k1,k2,iloc
      double precision :: xx,yy,zz,wx,wy,wz,wt(8),val1,val2,val3,val4
      logical :: notfindx,notfindy,notfindz
      iloc = 0
      notfindx =.true.
      do while(notfindx.and.i1.gt.0.and.i1.lt.nx)
         i2 = i1+1
         if(xx.lt.outx(i1))then
         i1 = i1-1
         else
            if(xx.le.outx(i2))then
            notfindx =.false.
            wx = (xx-outx(i1))/(outx(i2)-outx(i1))
            else
            i1 = i1+1
            end if
        end if
      end do
      notfindy =.true.
      do while(notfindy.and.j1.gt.0.and.j1.lt.ny)
         j2 = j1+1
         if(yy.lt.outy(j1))then
         j1 = j1-1
         else
            if(yy.le.outy(j2))then
            notfindy =.false.
            wy = (yy-outy(j1))/(outy(j2)-outy(j1))
            else
            j1 = j1+1
            end if
        end if
      end do
      notfindz =.true.
      do while(notfindz.and.k1.gt.0.and.k1.lt.nz)
         k2 = k1+1
         if(zz.lt.outz(k1))then
         k1 = k1-1
         else
            if(zz.le.outz(k2))then
            notfindz =.false.
            wz = (zz-outz(k1))/(outz(k2)-outz(k1))
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
      val1=wt(1)*vstart3d(i1,j1,k1)+wt(2)*vstart3d(i1,j2,k1)+wt(3)*vstart3d(i2,j2,k1)+&
           wt(4)*vstart3d(i2,j1,k1)+wt(5)*vstart3d(i1,j1,k2)+wt(6)*vstart3d(i1,j2,k2)+&
           wt(7)*vstart3d(i2,j2,k2)+wt(8)*vstart3d(i2,j1,k2)
      val2=wt(1)*vinvert3d(i1,j1,k1)+wt(2)*vinvert3d(i1,j2,k1)+wt(3)*vinvert3d(i2,j2,k1)+&
           wt(4)*vinvert3d(i2,j1,k1)+wt(5)*vinvert3d(i1,j1,k2)+wt(6)*vinvert3d(i1,j2,k2)+&
           wt(7)*vinvert3d(i2,j2,k2)+wt(8)*vinvert3d(i2,j1,k2)
      val3=wt(1)*outdv3d(i1,j1,k1)+wt(2)*outdv3d(i1,j2,k1)+wt(3)*outdv3d(i2,j2,k1)+&
           wt(4)*outdv3d(i2,j1,k1)+wt(5)*outdv3d(i1,j1,k2)+wt(6)*outdv3d(i1,j2,k2)+&
           wt(7)*outdv3d(i2,j2,k2)+wt(8)*outdv3d(i2,j1,k2)
      val4=wt(1)*outray3d(i1,j1,k1)+wt(2)*outray3d(i1,j2,k1)+wt(3)*outray3d(i2,j2,k1)+&
           wt(4)*outray3d(i2,j1,k1)+wt(5)*outray3d(i1,j1,k2)+wt(6)*outray3d(i1,j2,k2)+&
           wt(7)*outray3d(i2,j2,k2)+wt(8)*outray3d(i2,j1,k2)
      else
         print*,'cannot locate this point',xx,yy,zz,i1,j1,k1
         stop
      end if
      return
      end

  
      subroutine epidis(re,pe,rs,ps,dis)
      double precision :: pe,re,ps,rs,dis
      double precision :: pi,rad,p1,p2,r1,r2
      double precision :: xa,ya,za,xb,yb,zb,val
      dis = sqrt((re-rs)**2+(pe-ps)**2)
      return
      end



  


