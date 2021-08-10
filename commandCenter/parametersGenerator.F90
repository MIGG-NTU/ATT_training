
	program param
        implicit none
        include "../data/evtstaminmax"
	double precision,parameter :: eps = 1.0E-8
	double precision,parameter :: zbeg=-36.0,zend=2.0
	double precision,parameter :: dx=0.40,dy= 0.4,dz=0.4   ! x -> longitude, y -> latitude
	double precision :: xbeg,xend,ybeg,yend
	integer,parameter :: xnext=20, ynext=5
	integer :: nx,ny,nz

        ! (0) How many iterations will you conduct? (start from 1)
        integer,parameter :: niter = 1

        ! the maximum ratio of traveltime residual over obs and the absolute residual are bounded by tbd and tabs 
        double precision,parameter :: tbd=0.50,tabs=5.0

	double precision,parameter :: updatemax=0.02

        ! (6) how many sets of inversion grid for velocity inversion (you should prepare the meshes)
        integer,parameter :: ngrid = 5


        open(11,file='gridnumber')
        write(11,'(A26,I5)') "integer,parameter :: nset=",ngrid
        close(11)

        open(40,file='../inversion/parameter_inversion')
        write(40,'(A27,I6)') "integer,parameter :: niter=",niter
        write(40,'(A34,F12.5,A7,F12.5)') &
                "double precision,parameter :: tbd=",tbd,", tabs=",tabs 
        write(40,'(A40,F12.5)') "double precision,parameter :: updatemax=",updatemax
        close(40)


 	xbeg = xmin - xnext*dx
    	xend = xmax + xnext*dx
    	ybeg = ymin - ynext*dy
    	yend = ymax + ynext*dy    

   	nx = max(ifix(real((xend + eps - xbeg)/dx))+1,2)
  	ny = max(ifix(real((yend + eps - ybeg)/dy))+1,2)
   	nz = max(ifix(real((zend + eps - zbeg)/dz))+1,2)

        open(10,file='region') 
        write(10,'(A35,F12.5,A7,F12.5)') "double precision,parameter :: xbeg=",xbeg,", xend=",xbeg+(nx-1)*dx
        write(10,'(A35,F12.5,A7,F12.5)') "double precision,parameter :: ybeg=",ybeg,", yend=",ybeg+(ny-1)*dy
        write(10,'(A35,F12.5,A7,F12.5)') "double precision,parameter :: zbeg=",zbeg,", zend=",zbeg+(nz-1)*dz
        write(10,'(A33,F12.5,A5,F12.5,A5,F12.5)') "double precision,parameter :: dx=",dx,", dy=",dy,", dz=",dz
        write(10,'(A24,I5,A5,I5,A5,I5)') "integer,parameter :: nx=",nx,", ny=",ny,", nz=",nz
        close(10)


        stop
        end
