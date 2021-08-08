       module regmesh      
      integer :: ninv
      ! regular parameterization
      integer :: invx,invy,invz
      double precision,dimension(:),allocatable :: xinv,yinv,zinv
      double precision,dimension(:),allocatable :: att1d
      double precision,dimension(:,:),allocatable :: att2d
      double precision,dimension(:,:,:),allocatable :: att3d
      end module


      program meshgenerator
      integer :: meshtype
      call read_mesh_cube
      stop
      end


      subroutine read_mesh_cube
      use regmesh
      include '../commandCenter/region'
      include '../commandCenter/gridnumber'
      integer :: i,j,k,iset
      double precision :: ddx,ddy,ddz,disx,disy,disz
      character(len=3) :: cset
      character(len=12) :: filename
      disx = xend-xbeg
      disy = yend-ybeg
      disz = zend-zbeg
      invx = 14
      invy = 3 
      invz = 11 
      disx = disx/(invx-2)
      disy = disy/(invy-2)
      disz = disz/(invz-2)
      print*,'disx,disy,disz',disx,disy,disz
      ddx  = disx/nset
      ddy  = disy/nset
      ddz  = disz/nset
      allocate(xinv(invx))
      allocate(yinv(invy))
      allocate(zinv(invz))
      open(20,file='multiple-grid')
      write(20,*) nset,invx,invy,invz  
      do iset = 1,nset
         write(cset,'(I3.3)') iset
         filename = 'mesh_cube'//cset
         open(10,file=filename)
         do k=1,invz
            zinv(k)=zbeg+iset*ddz+(k-2)*disz
         end do
         do j=1,invy
            yinv(j)=ybeg+iset*ddy+(j-2)*disy
         end do
         do i=1,invx
            xinv(i)=xbeg+iset*ddx+(i-2)*disx
         end do
         xinv(1) =  xinv(1) - disx/2.0
      	xinv(invx) = xinv(invx) + disx/2.0
      	yinv(1) =  yinv(1) - disy/2.0
      	yinv(invy) = yinv(invy)  + disy/2.0
      	zinv(1) =  zinv(1) - 1.001
      	zinv(invz) = zinv(invz)  + 1.001
         yinv(1) = -10000.0
	      yinv(2) = 0
	      yinv(3) = 100000.0
      	write(10,*) invx,invy,invz
	      write(20,*) iset,invx,invy,invz
      	do k=1,invz
      	do j=1,invy
         do i=1,invx
            write(10,*) xinv(i),yinv(j),1.0*zinv(k)
            write(20,*) xinv(i),yinv(j),1.0*zinv(k)
      	end do
         end do
         end do
      	close(10)
      end do
      close(20)
      deallocate(xinv,yinv,zinv)
      return
      end
