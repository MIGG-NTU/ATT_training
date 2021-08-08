    program velocity3d_generate
    implicit none
    include '../commandCenter/region'
    integer :: ii, jj, kk
    double precision :: xx, yy, zz
    double precision,dimension(:,:,:),allocatable :: vel3d
 
    allocate(vel3d(nx,ny,nz))
    open(22,file='velocity3d')
    
    write(22,*) nx,ny,nz
    vel3d = 0
    do jj=1,ny
         yy   = ybeg+(jj-1)*dy
         do ii=1,nx
            xx   = xbeg+(ii-1)*dx 
            do kk=1,nz
               zz   = zbeg+(kk-1)*dz              
               vel3d(ii,jj,kk) = min(7.5,6.0-0.05*zz) 
               write(22,'(4F15.7)') xx,yy,zz,vel3d(ii,jj,kk)
            end do
        end do
    end do
    close(22)
    stop

    end
