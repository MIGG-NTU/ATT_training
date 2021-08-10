
	program param
	
        implicit none
        
        !*************************************  Part A: Select Your Research Area ******************************************
        !*
        !*              Please run workflowStep0. Once you are okay with the selected research area, please go to part B. 
        !*
        ! (1) Define your research region. Only check the earthquakes with magnitude between magmin and magmax. 
        !     If the magnitude greater than dumpmag, dump the earthquake as a big one.
        double precision,parameter :: xwest=0.0, xeast=240.0
        double precision,parameter :: ysouth=-1.0, ynorth=1.0
        double precision :: zbot=-40.0, ztop=5.0    ! z points upward
        double precision :: magmin=0.0,magmax=9.0,dumpmag=5.0
        double precision :: stwest=-1.0,steast=1.0,stsouth=-0.0,stnorth=0.0   
        !*********************  Part B: Select the traveltime data based on some criteria in your research area ****************
        !*
        ! (5) choose your interested phase 1 -> P , 2 -> S
        integer,parameter :: iphase=1
        ! (6) Choose the time window [idbeg,idend] YYMMDD
        integer,parameter :: idbeg = 19000101,idend=20220101
        !  size of searching block
        double precision,parameter :: dx=0.50, dy=0.50, ddep=0.50       
        integer,parameter :: NBD   = 1 ! The earthquake must have at least NBD phases to be selected. 

        open(12,file='earthquakeRegion')
	write(12,'(9F12.4)') xwest,xeast,ysouth,ynorth,zbot,ztop,magmin,magmax,dumpmag
        close(12)

 	open(13,file='stationRegion')
        write(13,*) xwest+stwest,xeast+steast,ysouth+stsouth,ynorth+stnorth
        close(13)


        open(14,file='phasewanted')
        write(14,'(A25,I10)') "integer,parameter :: ips=",iphase
        close(14)
      

        open(22,file='timePeriod')
        write(22,'(2I10)') idbeg,idend 
        close(22)

        open(23,file='blockSelection')
        write(23,'(3F12.5)') dx,dy,ddep
        write(23,'(I6)') nbd
        close(23)

        stop

        end

     
