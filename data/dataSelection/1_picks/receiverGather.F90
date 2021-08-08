	program receiverGathers
	include 'nevtstaray'
	integer :: irecsrc(nsta,nevt),irecord(nsta)
	double precision :: xsrc(nevt),ysrc(nevt),zsrc(nevt)
	double precision :: xrec(nsta),yrec(nsta),zrec(nsta)
	double precision :: tobs(nsta,nevt),tori(nevt)
	double precision :: ztmp,xtmp,ytmp,ttmp
	integer :: i,j,k

	tori=0.0

	open(10,file='receivers')
	do i=1,nsta
	   read(10,*) j,zrec(i),xrec(i),yrec(i)
	end do
	close(10)
	irecsrc = 0
	tobs    = 0
	irecord = 0
	open(20,file='obstime-rec-gather')
	do while(.true.)
	   read(20,*,end=999) i,j,k,ztmp,xtmp,ytmp,ttmp
	   xsrc(k)   = xtmp
           ysrc(k)   = ytmp
           zsrc(k)   = ztmp
	   tobs(i,k) = ttmp
	   irecsrc(i,k) = 1
	   irecord(i) = j
	end do
999	continue
	close(20)

	open(30,file='traveltimeReceiverGathers')
	do i=1,nsta
	   write(30, '(I9,3F16.5,I12)') i,zrec(i),xrec(i),yrec(i),irecord(i)
	   k = 0
	   do j=1,nevt
		if(irecsrc(i,j).gt.0)then
		   k = k + 1
		   write(30,'(3I9,4F16.5,F12.5)') i,k,j,zsrc(j),xsrc(j),ysrc(j),tori(j),tobs(i,j)
		end if 
	   end do
	end do
	
	close(30)	

	stop
	end
