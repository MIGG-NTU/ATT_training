    program velocity3d_true_generate
    implicit none
    include '../commandCenter/region'
    integer :: ii, jj, kk
    double precision :: xx, yy, zz, disx, disy, disz, dvoverv
    double precision,dimension(:,:,:),allocatable :: vel3d
 
    allocate(vel3d(nx,ny,nz))
    open(22,file='velocity3d_true')
    open(25,file='velocity2d_true')

    write(22,*) nx,ny,nz
    vel3d = 0
    do jj=1,ny
         yy   = ybeg+(jj-1)*dy
         disy = 2.0*acos(-1.0)*(1.0/100.0)*yy
         do ii=1,nx
            xx   = xbeg+(ii-1)*dx 
            disx = 2.0*acos(-1.0)*(1.0/80.0)*xx  
            do kk=1,nz
               zz   = zbeg+(kk-1)*dz
               if(zz.le.0.0)then
                  disz = 2.0*acos(-1.0)*(sqrt(49.0-8.0*zz)-7.0)/8
               else
                  disz = 0.0
               end if
               vel3d(ii,jj,kk) = min(7.5,6.0-0.05*zz) 
               dvoverv = 0.06*sin(disx)*sin(disz)
               write(22,'(4F15.7)') xx,yy,zz,vel3d(ii,jj,kk)*(1.0+dvoverv)
               if((yy.ge.-0.001).and.(yy.le.0.001))then
                  write(25,*) xx,zz,dvoverv
               end if
            end do
        end do
    end do
    close(22)
    close(25)
    stop

    end
