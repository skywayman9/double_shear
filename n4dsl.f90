!!!    This program solves viscous double periodic shear layer test by an adaptive upwind-central scheme

!!! 	Reconstruction schemes are well known 3rd order and 5th order upwind schemes (and the central 4th and 6th order schemes) but in a adaptive way.

!!! 	Componentwise Local Lax-Friedrich Riemann solver is used for inviscid fluxes

!!! 	Viscous fluxes are computed by the fourth order alpha-damping approach of Hiro Nishikawa. AIAA 2010-5093


!**********************************************************************
! Written by Sainath Ch, Email: sainath@caltech.edu
!**********************************************************************


	program shearlayer

	implicit none


    integer 				:: i, j
	integer					:: N, rk_step
	integer					:: NTMAX=100000

	double precision		:: t_end, time , dt

	double precision		:: CFL =0.4d0, gamma =(1.40d0)

	integer, parameter 		:: NX = 160, NY = 160, ghostp = 6, n_eqn = 4


	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, xmin, xmax, ymin, ymax

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn,0:2) ! 0:2 is for the three Runge-Kutta time steps


	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)

    double precision 		:: der_ux(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision 		:: der_vx(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision 		:: der_uy(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision 		:: der_vy(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    integer, parameter		:: file_save=600

	integer 				:: time_ini,time_end,time_calc
	double precision 		:: start, finish


	integer :: reconstruction


	common /grid/ dx, dy



	call cpu_time(start)
	write(*,*) 'Program start...'

	call system_clock(count = time_ini)

	write(*,*) "Enter 23 for 3rd/2nd order, 34 for 3rd/4th order or 56 for 5th/6th order reconstruction"

	read(*,*) reconstruction



		xmin = 0.0d0
		xmax = 1.0d0

		ymin = 0.0d0
		ymax = 1.0d0

		t_end = 1.00

	! Generate simple grid

		
	    dx = (xmax - xmin)/NX
		dy = (ymax - ymin)/NY

		do i = -ghostp, NX + ghostp

			  x(i) = xmin + (i-0.5d0)*dx
		enddo

		do j = -ghostp, NY + ghostp

		  	y(j) = ymin + (j-0.5d0)*dy
		
		enddo	

	call initialconditions(x,y,density, u_vel, v_vel,pressure, sound, gamma,NX,NY,ghostp)


	time = 0.0d0

	N=1
	
	! call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

	  write(*,*)'*********************************************'
      write(*,*)'   time step N        time             '
      write(*,*)'*********************************************'

    ! Computations starts here

	do while(time.lt.t_end)
			call timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)
			
			time = time + dt

			
			write(*,*) N ,time
			do rk_step = 0,2

				call boundaryconditions(density, u_vel, v_vel, pressure, x, y, time, gamma,NX,NY,ghostp)


    			do j = -ghostp, NY + ghostp
					do i = -ghostp, NX + ghostp

						cons(i,j,1,rk_step ) = density(i,j)
					    cons(i,j,2,rk_step ) = density(i,j)*u_vel(i,j)
					    cons(i,j,3,rk_step ) = density(i,j)*v_vel(i,j)
					    cons(i,j,4,rk_step ) = pressure(i,j)/(gamma-1.0d0) + 0.5d0*density(i,j)*(u_vel(i,j)**2.0d0 + v_vel(i,j)**2.0d0)

					enddo
				enddo

				call FX(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp,reconstruction)

				call GY(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp,reconstruction)

				call VF(density, u_vel, v_vel, pressure, residual, gamma, rk_step,NX,NY,ghostp)
				

				call rungekutta(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step, dt,NX,NY,ghostp,n_eqn)

			enddo

			N=N+1
			

			if(MOD(N,file_save) .eq. 0) then
				
				call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp,cons)
				
			endif	

			if (abs(time-t_end) .le. 1.0d-06) then
				
				call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp,cons)

				write(*,*)'*********************************************'
           		write(*,*)'   Number of time steps = ',N
          	    write(*,*)'*********************************************'

          	    exit
          	endif

    enddo
    	write(*,*)'*********************************************'
        write(*,*)'   Number of time steps = ',N
        write(*,*)'*********************************************'
    		

    	call boundaryconditions(density, u_vel, v_vel, pressure, x, y, time, gamma,NX,NY,ghostp)
    	call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp,cons)


	!***********************************************************************
	!*****                       For vorticity	                       *****
	!***********************************************************************

	   	do j=-3,NY+4
	        do i=-3,NX+4 
			der_vx(i,j)   =(-v_vel   (i-1,j) + v_vel   (i+1,j))/2.0d0
	        der_uy(i,j)   =(-u_vel   (i,j-1) + u_vel   (i,j+1))/2.0d0
	        enddo
	    enddo

    		call tecplot(N/file_save,density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp,der_vx,der_uy)

			call system_clock(count = time_end)
    
   			time_calc = time_end-time_ini
    
    	write(*,'(A20,I10,A)')   'Calculation time ',time_calc,' [CPU ticks]'

    	call cpu_time(finish)
    
    	write(*,*) " Total CPU time to solution = ", finish-start, " seconds"

    	write(*,*) 'Program ends...'


	end program shearlayer


	!***********************************************************************
	!*****                       Initial conditions                    *****
	!***********************************************************************


	subroutine initialconditions(x,y,density, u_vel, v_vel,pressure, sound, gamma,NX,NY,ghostp)

	implicit none

	integer 				:: i, j

	integer			 		:: NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, xmin, xmax, ymin, ymax, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision, parameter :: Pi = ACOS(-1.d0), delW = 120.d0, delP = 0.05d0


	do j = 1, NY
		do i = 1, NX
		
			if (y(j) .le. 0.5d0) then
				u_vel(i,j) = TANH(delW*(y(j)-0.25d0))
			else
				u_vel(i,j) = TANH(delW*(0.75d0-y(j)))
			endif
		
			v_vel(i,j) 		= delP*SIN(2.d0*Pi*(x(i)+0.0d0))
			pressure(i,j) 	= 1.d0/(gamma*0.1d0**2)
			density(i,j)    =  gamma*0.1d0**2*pressure(i,j)
			
		enddo
	enddo
	

	end subroutine initialconditions

	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp,cons)

	integer 				:: i, j
	integer			 		:: NX, NY, ghostp
	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4,0:2)

	open(unit=25, file="solnsd.txt",action="write",status="replace")
	do j = 1, NY
	    do i = 1, NX
	 
	  	write(25,'(7F25.8)') x(i),y(j),density(i,j),pressure(i,j),u_vel(i,j),v_vel(i,j),cons(i,j,4,2)
	 	
	 	enddo
	enddo
	close(25)


 	end subroutine output

	!***********************************************************************
	!*****                      TEC Output 			                   *****
	!***********************************************************************

 	subroutine tecplot(file,density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp,der_vx,der_uy)
	implicit none


    integer 				:: i, j, file,l
	integer			 		:: NX, NY, ghostp
	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	character(len=8) 		:: number*4, file_name

	double precision 		:: der_vx(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision 		:: der_uy(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

		write(number,'(i4.4)') file
		file_name="Rslt"//number
		open(unit=1,file=file_name//'.plt')
		
		write(1,*) 'TITLE="',file_name,'"'
	    write(1,*) 'VARIABLES = "x","y","rho","vx","vy","Pre"'
		write(1,*) "ZONE I=",NX," J=",NY," F=POINT"


		    do j = 1, NY
		        do i = 1, NX

		          write(1,'(6F25.8)') x(I), y(J), density(I,J), u_vel(I,J), v_vel(I,J), (der_vx(i,j)-der_uy(i,j))

			  	enddo
		    enddo
		
	    close(1)


    end subroutine tecplot


 	!***********************************************************************
	!*****                       Compute time step                	   *****
	!***********************************************************************

 	subroutine timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)
 	implicit none


 	integer 				:: i, j

	integer			 		:: NX, NY, ghostp
	double precision 		:: dx, dy,gamma

	double precision 		:: dt, time, t_end, dtnew, CFL
 	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
 	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: x_velocity, y_velocity

	double precision        ::mu_lam,dt_visc
	common /grid/ dx, dy

	mu_lam  = 1.0d0/10000.0d0
	dt_visc = 0.375d0*DMIN1( dx**2/mu_lam,  dy**2/mu_lam)

	dt = 1.0d10

		do i = 1, NX

			do j = 1, NY

				sound(i,j) =  (gamma*pressure(i,j)/density(i,j))**0.5d0
				x_velocity =  ABS(u_vel(i,j)) + sound(i,j)
				y_velocity =  ABS(v_vel(i,j)) + sound(i,j)

				dtnew = min(dx/x_velocity, dy/y_velocity)

				if(dtnew .lt. dt) dt = DMIN1(dtnew, dt_visc)
			
			enddo
		
		enddo

		dt = CFL*dt

		if ((time+dt) .gt. t_end ) then

			dt = t_end - time
		
		endif	


 	end subroutine timestep


 	!***********************************************************************
	!*****                       Boundary conditions - Peri            *****
	!***********************************************************************


 	subroutine boundaryconditions(density, u_vel, v_vel, pressure, x, y, time, gamma,NX,NY,ghostp)
    implicit none

	integer 				:: i, j, NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, x1, y1, gamma

	double precision		:: time

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision, parameter 		:: pi=acos(-1.0d0),period = 30.d0/2.68d0

	common /mesh/ dx, dy

		! Set left and right boundary condition

			do i = 1,ghostp
				
				do j = 1, NY

					density(-i+1,j)	    = density(NX-i+1,j)
					u_vel(-i+1,j)		= u_vel (NX-i+1,j)
					v_vel(-i+1,j)		= v_vel (NX-i+1,j)
					pressure(-i+1,j)	= pressure(NX-i+1,j)


					density(NX+i,j)		= density(i,j)
					u_vel(NX+i,j)		= u_vel(i,j)
					v_vel(NX+i,j)		= v_vel(i,j)
					pressure(NX+i,j)	= pressure(i,j)

				enddo
			
			enddo

	! Set top and bottom boundary condition

			do j = 1,ghostp
				
				do i = 1, NX

					density(i,-j+1)  = density(i,NY-j+1)
					u_vel(i,-j+1)    = u_vel(i,NY-j+1)
					v_vel(i,-j+1)    = v_vel(i,NY-j+1)
					pressure(i,-j+1) = pressure(i,NY-j+1)


					density(i,NY+j)		=  density(i,j)
					u_vel(i,NY+j)		=  u_vel(i,j)
					v_vel(i,NY+j)		=  v_vel(i,j)
					pressure(i,NY+j)	=  pressure(i,j)

				enddo
			
			enddo
    
    
     
      end subroutine boundaryconditions

 	!***********************************************************************
	!*****                       Time step, TVD- Runge Kutta           *****
	!***********************************************************************

	subroutine rungekutta(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step, dt,NX,NY,ghostp,n_eqn)
	implicit none

	integer 				:: i, j, k, rk_step, rk

	integer 				:: NX, NY, ghostp,n_eqn

	double precision 		:: dt, uv, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn,0:2) ! 0:2 is for the three Runge-Kutta time steps

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)


		if(rk_step .EQ.0) then

			do k=1, n_eqn

			    do j = 1, NY
			
				  do i = 1, NX
			
				    cons(i,j,k,1) = cons(i,j,k,0) + dt*residual(i,j,k)
			
				  enddo
			
				enddo
			enddo

		elseif(rk_step .EQ.1) then

			do k=1, n_eqn
			    
			    do j = 1, NY
				
				  do i = 1, NX
				
				    cons(i,j,k,2) = (3.0d0/4.0d0)*cons(i,j,k,0) + (1.0/4.0d0)*(cons(i,j,k,1) + dt*residual(i,j,k))
				
				  enddo
				
				enddo
			enddo

		else

			do k=1, n_eqn
			    
			    do j = 1, NY
				
				  do i = 1, NX
				
				    cons(i,j,k,0) = (1.0d0/3.0d0)*cons(i,j,k,0) + (2.0d0/3.0d0)*(cons(i,j,k,2) + dt*residual(i,j,k))
				
				  enddo
				
				enddo
			
			enddo

		endif
			
				rk = MOD(rk_step +1, 3)

					do i = 1, NX
						do j = 1, NY

					    density(i,j)		= cons(i,j,1,rk)
				       	u_vel(i,j)			= cons(i,j,2,rk)/density(i,j)
				        v_vel(i,j)			= cons(i,j,3,rk)/density(i,j)
				        uv 					= u_vel(i,j)**2 + v_vel(i,j)**2
					    pressure(i,j)		= (gamma-1.0)*(cons(i,j,4,rk)-0.5*cons(i,j,1,rk)*uv)


						enddo
					enddo

	end subroutine rungekutta




 	subroutine FX(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp,reconstruction)
 	implicit none


 	integer     		:: i, j, rk_step, M, NX, NY, ghostp,k
 	integer, parameter	::  n_eqn =4

 	double precision    ::	 dx, dy, gamma

 	integer,intent(IN)  :: reconstruction

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: sound(-ghostp:NX+ghostp, -ghostp:NY+ghostp),enthalpy(-ghostp:NX+ghostp)
	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4,0:2),constime(-ghostp:NX+ghostp,4), primitive(-ghostp:NX+ghostp,4) ! 0:2 is for the three Runge-Kutta time steps

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4)

	double precision		:: lefteigen(-ghostp:NX+ghostp,4,4) , righteigen(-ghostp:NX+ghostp,4,4) 		! 2D 4 by 4 matrix 4 eigen values and .....
	double precision		:: consl(-ghostp:NX+ghostp,4),consr(-ghostp:NX+ghostp,4)
	double precision		:: priml(-ghostp:NX+ghostp,4),primr(-ghostp:NX+ghostp,4)

	double precision		:: flux_half(-ghostp:NX+ghostp,4)
	double precision		:: fright(-ghostp:NX+ghostp,4), fleft(-ghostp:NX+ghostp,4)

	double precision		:: lambda1(-ghostp:NX+ghostp),lambda2(-ghostp:NX+ghostp),lambda3(-ghostp:NX+ghostp),lambda4(-ghostp:NX+ghostp)
	double precision		:: den_average(-ghostp:NX+ghostp)
	double precision   		:: enthalpy_average, uv_average, sound_average, u_average, v_average,sound_left,sound_right
	double precision   		:: T0, T1, T2, T3, B1, B2,inverse_sound,Da, ML,MR

	!!!! Riemann solver

	double precision		:: densityleft(-ghostp:NX+ghostp), pressureleft(-ghostp:NX+ghostp),u_velleft(-ghostp:NX+ghostp), v_velleft(-ghostp:NX+ghostp)
	double precision		:: densityright(-ghostp:NX+ghostp), pressureright(-ghostp:NX+ghostp),u_velright(-ghostp:NX+ghostp), v_velright(-ghostp:NX+ghostp)
	double precision		:: enthalpyleft(-ghostp:NX+ghostp), enthalpyright(-ghostp:NX+ghostp),energyleft(-ghostp:NX+ghostp),energyright(-ghostp:NX+ghostp)

	double precision 		:: delrho(-ghostp:NX+ghostp), delmu(-ghostp:NX+ghostp),delmv(-ghostp:NX+ghostp)
	double precision  		:: deltoten(-ghostp:NX+ghostp)
	double precision 		:: dissipation(-ghostp:NX+ghostp,4)


	double precision		:: alpha1,alpha2,alpha3,alpha4
	double precision 		:: SL,SR,SP,EL,ER
	double precision		:: mx,my

	double precision		:: sqrt_rho_R,sqrt_rho_L,rho,divisor
    double precision		:: uv_left(-ghostp:NX+ghostp),uv_right(-ghostp:NX+ghostp)


	common /grid/ dx, dy



	mx=1.0d0;my=0.0d0
	

	do j = 1, NY

		do i = -ghostp, NX + ghostp


		    constime(i,1) = density(i,j)
		    constime(i,2) = density(i,j)*u_vel(i,j)
		    constime(i,3) = density(i,j)*v_vel(i,j)
		    constime(i,4) = pressure(i,j)/(gamma-1.0d0) + 0.5d0*density(i,j)*(u_vel(i,j)**2.0d0 + v_vel(i,j)**2.0d0)

		enddo
	  
	  
	  call recon(NX, constime, consl,consr,mx,my,ghostp,n_eqn,reconstruction)


	  	do i =-1,NX+1


			densityleft(i)     = consl(i,1)
			u_velleft(i)  	   = consl(i,2)/densityleft(i)
			v_velleft(i)	   = consl(i,3)/densityleft(i)

			densityright(i)     = consr(i,1)
			u_velright(i)   	= consr(i,2)/densityright(i)
			v_velright (i)      = consr(i,3)/densityright(i)




			densityleft(i)     = consl(i,1)
			u_velleft(i)  	   = consl(i,2)/densityleft(i)
			v_velleft(i)	   = consl(i,3)/densityleft(i)
			uv_left(i)		   = u_velleft(i)**2 + v_velleft(i)**2
			pressureleft (i)   = (gamma-1.0)*(consl(i,4)-0.5*consl(i,1)*uv_left(i))
			enthalpyleft(i)	   = (pressureleft(i) + consl(i,4))/densityleft(i)

			sound_left		   = (gamma*pressureleft(i)/densityleft(i))**0.5d0


			densityright(i)     = consr(i,1)
			u_velright(i)   	= consr(i,2)/densityright(i)
			v_velright (i)      = consr(i,3)/densityright(i)
			uv_right(i)		    = u_velright(i)**2 + v_velright(i)**2
			pressureright(i)    = (gamma-1.0)*(consr(i,4)-0.5*consr(i,1)*uv_right(i))
			enthalpyright(i)    = (pressureright(i) + consr(i,4))/densityright(i) 

			sound_right		   = (gamma*pressureright(i)/densityright(i))**0.5d0



	  		fleft(i,1)	   = consl(i,2)
			fleft(i,2)	   = consl(i,2) * u_velleft(i) + pressureleft(i)
			fleft(i,3)	   = consl(i,2) * v_velleft(i) 
			fleft(i,4)	   = u_velleft(i) * (consl(i,4) + pressureleft(i))

			fright(i,1)	   = consr(i,2)
			fright(i,2)	   = consr(i,2) * u_velright(i)+ pressureright(i)
			fright(i,3)	   = consr(i,2) * v_velright(i) 
			fright(i,4)	   = u_velright(i) * (consr(i,4) + pressureright(i))

			delrho(i)			= consr(i,1)-consl(i,1)
			delmu(i)			= consr(i,2)-consl(i,2)
			delmv(i)			= consr(i,3)-consl(i,3)
			deltoten(i)			= consr(i,4)-consl(i,4)



				sqrt_rho_L = sqrt(densityleft(i))
				sqrt_rho_R = sqrt(densityright(i))

				rho 	   = sqrt(densityright(i)/densityleft(i))*densityleft(i)

				divisor	   = 1.0d0/(sqrt_rho_R+sqrt_rho_L)

				u_average 		= (  (u_velleft(i)*sqrt_rho_L)    + (   u_velright(i)*sqrt_rho_R))*divisor
				v_average 		= (  (v_velleft(i)*sqrt_rho_L)    + (   v_velright(i)*sqrt_rho_R))*divisor
				enthalpy_average= (  (enthalpyleft(i)*sqrt_rho_L) + (enthalpyright(i)*sqrt_rho_R))*divisor
				uv_average 		= 	0.5d0 * (u_average**2.0d0 + v_average**2.0d0)
				sound_average 	= 	sqrt((gamma - 1.0d0) * (enthalpy_average - uv_average))

				T0 					= 	u_average * sound_average


				righteigen(i,1,1) = 1.0d0
				righteigen(i,1,2) = 0.0d0
				righteigen(i,1,3) = 1.0d0
				righteigen(i,1,4) = 1.0d0

				righteigen(i,2,1) = u_average - sound_average
				righteigen(i,2,2) = 0.0
				righteigen(i,2,3) = u_average
				righteigen(i,2,4) = u_average + sound_average

				righteigen(i,3,1) = v_average
				righteigen(i,3,2) = 1.0
				righteigen(i,3,3) = v_average
				righteigen(i,3,4) = v_average

				righteigen(i,4,1) = enthalpy_average - T0
				righteigen(i,4,2) = v_average
				righteigen(i,4,3) = uv_average
				righteigen(i,4,4) = enthalpy_average + T0

				inverse_sound 	= 1.0/sound_average
				B1 				= (gamma - 1.0) * inverse_sound**2
				B2 				= uv_average * B1
				T0 				= u_average * inverse_sound
				T1 				= B1 * u_average
				T2 				= 0.5 * B1
				T3 				= B1 * v_average


				lefteigen(i,1,1) = 0.5 * (B2 + T0)
				lefteigen(i,1,2) = - 0.5 * (T1 + inverse_sound)
				lefteigen(i,1,3) = - 0.5 * T3
				lefteigen(i,1,4) = T2

				lefteigen(i,2,1) = - v_average
				lefteigen(i,2,2) = 0.0
				lefteigen(i,2,3) = 1.0
				lefteigen(i,2,4) = 0.0

				lefteigen(i,3,1) = 1.0 - B2
				lefteigen(i,3,2) = T1
				lefteigen(i,3,3) = T3
				lefteigen(i,3,4) = - B1

				lefteigen(i,4,1) = 0.5 * (B2 - T0)
				lefteigen(i,4,2) = - 0.5 * (T1 - inverse_sound)
				lefteigen(i,4,3) = - 0.5 * T3
				lefteigen(i,4,4) = T2

				alpha1 			=  lefteigen(i,1,1)*delrho(i) + lefteigen(i,1,2)*delmu(i) + &
									   lefteigen(i,1,3)*delmv(i)  + lefteigen(i,1,4)*deltoten(i) 

				alpha2			= lefteigen(i,2,1)*delrho(i)  + lefteigen(i,2,2)*delmu(i) +&
								  lefteigen(i,2,3)*delmv(i)   + lefteigen(i,2,4)*deltoten(i)

				alpha3 			= lefteigen(i,3,1)*delrho(i) + lefteigen(i,3,2)*delmu(i) + &
								  lefteigen(i,3,3)* delmv(i) + lefteigen(i,3,4)*deltoten(i)

				alpha4			= lefteigen(i,4,1)*delrho(i) + lefteigen(i,4,2)*delmu(i) + &
								  lefteigen(i,4,3)* delmv(i) + lefteigen(i,4,4)*deltoten(i)



				!Componentwise LLF

				! Fleischmann, Nico, et al. "A low dissipation method to cure the grid-aligned shock instability." Journal of Computational Physics 401 (2020): 109004.

				lambda1(i) =  max(abs(u_velleft(i)-sound_left),abs(u_velright(i)-sound_right))
				lambda2(i) =  max(abs(u_velleft(i)),abs(u_velright(i)))
				lambda3(i) =  lambda2(i)
				lambda4(i) =  max(abs(u_velleft(i)+sound_left),abs(u_velright(i)+sound_right))



				dissipation(i,1) = (lambda1(i)*alpha1*righteigen(i,1,1) + lambda2(i)*alpha2*righteigen(i,1,2) &
									 	+ lambda3(i)*alpha3*righteigen(i,1,3)+ lambda4(i)*alpha4*righteigen(i,1,4))

				dissipation(i,2) = (lambda1(i)*alpha1*righteigen(i,2,1) + lambda2(i)*alpha2*righteigen(i,2,2) &
				                   		+ lambda3(i)*alpha3*righteigen(i,2,3)+ lambda4(i)*alpha4*righteigen(i,2,4))

				dissipation(i,3) =	((lambda1(i)*alpha1*righteigen(i,3,1) + &
									   		lambda2(i)*alpha2*righteigen(i,3,2) + &
									   		lambda3(i)*alpha3*righteigen(i,3,3))) + lambda4(i)*alpha4*righteigen(i,3,4)

				dissipation(i,4) =	((lambda1(i)*alpha1*righteigen(i,4,1) + &
									   		lambda2(i)*alpha2*righteigen(i,4,2) + &
									   		lambda3(i)*alpha3*righteigen(i,4,3))) + lambda4(i)*alpha4*righteigen(i,4,4)



				flux_half(i,:) =  0.5d0*((fright(i,:) +fleft(i,:) - dissipation(i,:)))	




		enddo


	    do i = 1, NX
		  residual(i,j,:) = -(flux_half(i,:) - flux_half(i-1,:))/dx
		enddo
	
	enddo
 	

 	end subroutine FX

 	!***********************************************************************
	!*****                      Flux in Y-direction                    *****
	!***********************************************************************

	subroutine GY(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step,NX,NY,ghostp,reconstruction)
	implicit none

 	integer     		:: i, j, rk_step, M, NX, NY, ghostp,k
 	integer, parameter	::  n_eqn =4

 	integer,intent(IN) :: reconstruction

 	double precision    ::	 dx, dy, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: sound(-ghostp:NX+ghostp, -ghostp:NY+ghostp),enthalpy(-ghostp:NY+ghostp)
	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4,0:2),constime(-ghostp:NY+ghostp,4),primitive(-ghostp:NY+ghostp,4) ! 0:2 is for the three Runge-Kutta time steps
	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4)

	double precision		:: lefteigen(-ghostp:NY+ghostp,4,4) , righteigen(-ghostp:NY+ghostp,4,4) 		! 2D 4 by 4 matrix 4 eigen values and .....

	double precision		:: consl(-ghostp:NY+ghostp,4),consr(-ghostp:NY+ghostp,4)
	double precision		:: priml(-ghostp:NY+ghostp,4),primr(-ghostp:NY+ghostp,4)

	double precision		:: flux_half(-ghostp:NY+ghostp,4)
	double precision		:: fright(-ghostp:NY+ghostp,4), fleft(-ghostp:NY+ghostp,4)

	double precision		:: lambda1(-ghostp:NY+ghostp),lambda2(-ghostp:NY+ghostp),lambda3(-ghostp:NY+ghostp),lambda4(-ghostp:NY+ghostp)
	double precision		:: den_average(-ghostp:NY+ghostp)
	double precision   		:: enthalpy_average, uv_average, sound_average, u_average, v_average, sound_right,sound_left
	double precision   		:: T0, T1, T2, T3, B1, B2,inverse_sound,Da,ML,MR

	!!!! Riemann solver

	double precision		:: densityleft(-ghostp:NY+ghostp), pressureleft(-ghostp:NY+ghostp),u_velleft(-ghostp:NY+ghostp), v_velleft(-ghostp:NY+ghostp)
	double precision		:: densityright(-ghostp:NY+ghostp), pressureright(-ghostp:NY+ghostp),u_velright(-ghostp:NY+ghostp), v_velright(-ghostp:NY+ghostp)
	double precision		:: enthalpyleft(-ghostp:NY+ghostp), enthalpyright(-ghostp:NY+ghostp),energyleft(-ghostp:NY+ghostp),energyright(-ghostp:NY+ghostp)


	double precision 		:: delrho(-ghostp:NY+ghostp), delmu(-ghostp:NY+ghostp),delmv(-ghostp:NY+ghostp)
	double precision  		:: deltoten(-ghostp:NY+ghostp)
	double precision 		:: dissipation(-ghostp:NY+ghostp,4)

	double precision		:: alpha1,alpha2,alpha3,alpha4
	double precision 		:: SL,SR,SP,EL,ER
	double precision		:: mx,my
	
	double precision		:: sqrt_rho_R,sqrt_rho_L,rho,divisor
	double precision		:: uv_left(-ghostp:NY+ghostp),uv_right(-ghostp:NY+ghostp)




	common /grid/ dx, dy


	mx=0.0d0;my=1.0d0

	

	do j = 1, NX

	 		do i = -ghostp, NY + ghostp


			   	constime(i,1) = density(j,i)
			    constime(i,2) = density(j,i)*u_vel(j,i)
			    constime(i,3) = density(j,i)*v_vel(j,i)
			    constime(i,4) = pressure(j,i)/(gamma-1.0d0) + 0.5d0*density(j,i)*(u_vel(j,i)**2.0d0 + v_vel(j,i)**2.0d0)

			enddo	  


	  call recon(NY, constime, consl,consr,mx,my,ghostp,n_eqn,reconstruction)

! Calculations for the Riemann solver

	 	do i =-1, NY+1

			densityleft(i)     = consl(i,1)
			u_velleft(i)  	   = consl(i,2)/densityleft(i)
			v_velleft(i)	   = consl(i,3)/densityleft(i)


			densityright(i)     = consr(i,1)
			u_velright(i)   	= consr(i,2)/densityright(i)
			v_velright (i)      = consr(i,3)/densityright(i)


			densityleft(i)     = consl(i,1)
			u_velleft(i)  	   = consl(i,2)/densityleft(i)
			v_velleft(i)	   = consl(i,3)/densityleft(i)
			uv_left(i)		   = u_velleft(i)**2 + v_velleft(i)**2
			pressureleft (i)   = (gamma-1.0)*(consl(i,4)-0.5*consl(i,1)*uv_left(i))
			enthalpyleft(i)	   = (pressureleft(i) + consl(i,4))/densityleft(i)

			sound_left		   = (gamma*pressureleft(i)/densityleft(i))**0.5d0



			densityright(i)     = consr(i,1)
			u_velright(i)   	= consr(i,2)/densityright(i)
			v_velright (i)      = consr(i,3)/densityright(i)
			uv_right(i)		    = u_velright(i)**2 + v_velright(i)**2
			pressureright(i)    = (gamma-1.0)*(consr(i,4)-0.5*consr(i,1)*uv_right(i))
			enthalpyright(i)    = (pressureright(i) + consr(i,4))/densityright(i) 

			sound_right		   = (gamma*pressureright(i)/densityright(i))**0.5d0



			fleft(i,1)	   = consl(i,3)
			fleft(i,2)	   = consl(i,2) * v_velleft(i)
			fleft(i,3)	   = consl(i,3) * v_velleft(i) + pressureleft(i)
			fleft(i,4)	   = consl(i,3)*enthalpyleft(i)

			fright(i,1)	   = consr(i,3)
			fright(i,2)	   = consr(i,2) * v_velright(i)
			fright(i,3)	   = consr(i,3) * v_velright(i) + pressureright(i)
			fright(i,4)	   = consr(i,3)*enthalpyright(i)



			delrho(i)			= consr(i,1)-consl(i,1)
			delmu(i)			= consr(i,2)-consl(i,2)
			delmv(i)			= consr(i,3)-consl(i,3)
			deltoten(i)			= consr(i,4)-consl(i,4)


			   	sqrt_rho_L = sqrt(densityleft(i))
				sqrt_rho_R = sqrt(densityright(i))

				rho 	   = sqrt(densityright(i)/densityleft(i))*densityleft(i)

				divisor	   = 1.0d0/(sqrt_rho_R+sqrt_rho_L)

				u_average 		= (  (u_velleft(i)*sqrt_rho_L)    + (   u_velright(i)*sqrt_rho_R))*divisor
				v_average 		= (  (v_velleft(i)*sqrt_rho_L)    + (   v_velright(i)*sqrt_rho_R))*divisor
				enthalpy_average= (  (enthalpyleft(i)*sqrt_rho_L) + (enthalpyright(i)*sqrt_rho_R))*divisor
				uv_average 		= 	0.5d0 * (u_average**2.0d0 + v_average**2.0d0)
				sound_average 	= 	sqrt((gamma - 1.0d0) * (enthalpy_average - uv_average))		
			    T0 				= 	v_average * sound_average

			   	righteigen(i,1,1) = 1.0
			    righteigen(i,1,2) = 0.0
			    righteigen(i,1,3) = 1.0
			    righteigen(i,1,4) = 1.0

			    righteigen(i,2,1) = u_average
			    righteigen(i,2,2) = 1.0
			    righteigen(i,2,3) = u_average
			    righteigen(i,2,4) = u_average

			    righteigen(i,3,1) = v_average - sound_average
			    righteigen(i,3,2) = 0.0
			    righteigen(i,3,3) = v_average
			    righteigen(i,3,4) = v_average + sound_average

			    righteigen(i,4,1) = enthalpy_average - T0
			    righteigen(i,4,2) = u_average
			    righteigen(i,4,3) = uv_average
			    righteigen(i,4,4) = enthalpy_average + T0

			    inverse_sound 		= 1.0/sound_average
			    B1 					= (gamma - 1.0) * inverse_sound**2
			    B2 					= uv_average * B1
			    T0 					= v_average * inverse_sound
			    T1 					= B1 * v_average
			    T2 					= 0.5 * B1
			    T3 					= B1 * u_average

			    lefteigen(i,1,1) = 0.5 * (B2 + T0)
			    lefteigen(i,1,2) = - 0.5 * T3
			    lefteigen(i,1,3) = - 0.5 * (T1 + inverse_sound)
			    lefteigen(i,1,4) = T2

			    lefteigen(i,2,1) = - u_average
			    lefteigen(i,2,2) = 1.0
			    lefteigen(i,2,3) = 0.0
			    lefteigen(i,2,4) = 0.0

			    lefteigen(i,3,1) = 1.0 - B2
			    lefteigen(i,3,2) = T3
			    lefteigen(i,3,3) = T1
			    lefteigen(i,3,4) = - B1

			    lefteigen(i,4,1) = 0.5 * (B2 - T0)
			    lefteigen(i,4,2) = - 0.5 * T3
			    lefteigen(i,4,3) = - 0.5 * (T1 - inverse_sound)
			    lefteigen(i,4,4) = T2


			        alpha1 			=  lefteigen(i,1,1)*delrho(i) + lefteigen(i,1,2)*delmu(i) + &
		     						   lefteigen(i,1,3)*delmv(i)  + lefteigen(i,1,4)*deltoten(i) 

		     		alpha2			= lefteigen(i,2,1)*delrho(i)  + lefteigen(i,2,2)*delmu(i) +&
		          					  lefteigen(i,2,3)*delmv(i)   + lefteigen(i,2,4)*deltoten(i)

		          	alpha3 			= lefteigen(i,3,1)*delrho(i) + lefteigen(i,3,2)*delmu(i) + &
		           					  lefteigen(i,3,3)* delmv(i) + lefteigen(i,3,4)*deltoten(i)

		           	alpha4			= lefteigen(i,4,1)*delrho(i) + lefteigen(i,4,2)*delmu(i) + &
		           					  lefteigen(i,4,3)* delmv(i) + lefteigen(i,4,4)*deltoten(i)

		!Componentwise LLF

		! Fleischmann, Nico, et al. "A low dissipation method to cure the grid-aligned shock instability." Journal of Computational Physics 401 (2020): 109004.
							

				    lambda1(i) =  max(abs(v_velleft(i)-sound_left),abs(v_velright(i)-sound_right))
					lambda2(i) =  max(abs(v_velleft(i)),abs(v_velright(i)))
					lambda3(i) =  lambda2(i)
					lambda4(i) =  max(abs(v_velleft(i)+sound_left),abs(v_velright(i)+sound_right))


				
				dissipation(i,1) = (lambda1(i)*alpha1*righteigen(i,1,1) + lambda2(i)*alpha2*righteigen(i,1,2) &
										 	+ lambda3(i)*alpha3*righteigen(i,1,3)+ lambda4(i)*alpha4*righteigen(i,1,4))

				dissipation(i,2) = (lambda1(i)*alpha1*righteigen(i,2,1) + lambda2(i)*alpha2*righteigen(i,2,2) &
					                   		+ lambda3(i)*alpha3*righteigen(i,2,3)+ lambda4(i)*alpha4*righteigen(i,2,4))
				
				dissipation(i,3) =	((lambda1(i)*alpha1*righteigen(i,3,1) + &
		     						   		lambda2(i)*alpha2*righteigen(i,3,2) + &
		     						   		lambda3(i)*alpha3*righteigen(i,3,3))) + lambda4(i)*alpha4*righteigen(i,3,4)

				dissipation(i,4) =	((lambda1(i)*alpha1*righteigen(i,4,1) + &
		     						   		lambda2(i)*alpha2*righteigen(i,4,2) + &
		     						   		lambda3(i)*alpha3*righteigen(i,4,3))) + lambda4(i)*alpha4*righteigen(i,4,4)



				flux_half(i,:) =  0.5d0*((fright(i,:) +fleft(i,:) - dissipation(i,:)))


		enddo


	    do i = 1, NY
	      residual(j,i,:) = residual(j,i,:) - (flux_half(i,:) - flux_half(i-1,:))/dy
		enddo

	enddo


	end subroutine GY

	!***********************************************************************
	!*****                      Viscous Flux in X-Y-directions         *****
	!***********************************************************************


	! Nishikawa, Hiroaki. "Beyond interface gradient: a general principle for constructing diffusion schemes." 40th fluid dynamics conference and exhibit. 2010.

    subroutine	VF(density, u_vel, v_vel, pressure, residual, gamma, rk_step,NX,NY,ghostp)
	
	implicit none
	

	integer :: i,j,k, NX, NY, ghostp,rk_step
 	
 	integer, parameter		::  n_eqn =4

 	double precision    	::	 dx, dy
	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp),temperature(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: mu_lam,kappa,gamma,Rgas,Pr


	double precision :: fluxvX(0:NX,0:NY,4),fluxvY(0:NX,0:NY,4)
	double precision :: derua(2), derva(2), derwa(2), derTa(2)
	double precision :: delV, vflux(4), sf(3)

    double precision :: der_ux(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: der_vx(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: der_uy(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: der_vy(-ghostp:NX+ghostp,-ghostp:NY+ghostp)


	double precision :: derT(-ghostp:NX+ghostp,-ghostp:NY+ghostp,1:2)
	

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)
	double precision 		:: u_left, u_right, v_left,v_right, p_left, p_right, d_left, d_right, p_avg, d_avg


    double precision :: vfx(4)
    double precision :: Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz,Theta(3), uface(3)



	common /grid/ dx, dy
   

    
    Pr= 0.73d0


	mu_lam = 1.0d0/10000.0d0
	kappa = gamma*mu_lam/(Pr*0.4d0)
	Rgas = 1.d0

    do j=-3,NY+4
        do i=-3,NX+4 
		
	        der_ux(i,j)   =(-u_vel   (i-1,j) + u_vel   (i+1,j))/2.0d0
	        der_vx(i,j)   =(-v_vel   (i-1,j) + v_vel   (i+1,j))/2.0d0

	        der_uy(i,j)   =(-u_vel   (i,j-1) + u_vel   (i,j+1))/2.0d0
	        der_vy(i,j)   =(-v_vel   (i,j-1) + v_vel   (i,j+1))/2.0d0

        enddo
    enddo

	do j=-3,NY+4
   		do i=-3,NX+4
        	temperature(i,j) = pressure(i,j) / (density(i,j)*Rgas)
    	enddo
    enddo

    do j=-1,NY+2
        do i=-1,NX+2 
		
        derT(i,j,1)   =(-temperature   (i-1,j) + temperature   (i+1,j))/2.0d0

        derT(i,j,2)   =(-temperature   (i,j-1) + temperature   (i,j+1))/2.0d0
        enddo
    enddo



! flux on I faces
    sf(1)=1; sf(2)=0
     do j = 1, NY
    	do i = 0, NX


    	u_left  = u_vel(i+0,j)+(1.0d0/2.0d0)*der_ux(i+0,j)
        u_right = u_vel(i+1,j)-(1.0d0/2.0d0)*der_ux(i+1,j)

        v_left  = v_vel(i+0,j)+(1.0d0/2.0d0)*der_vx(i+0,j)
        v_right = v_vel(i+1,j)-(1.0d0/2.0d0)*der_vx(i+1,j)

		derua(1) = (0.5d0*(der_ux(i,j) + der_ux(i+1,j))+(4.0d0/(3.0d0))*(u_right-u_left))/dx
		derva(1) = (0.5d0*(der_vx(i,j) + der_vx(i+1,j))+(4.0d0/(3.0d0))*(v_right-v_left))/dx

		derua(2) = (0.5d0*(der_uy(i,j) + der_uy(i+1,j))+(4.0d0/(3.0d0))*(u_right-u_left))/dx
		derva(2) = (0.5d0*(der_vy(i,j) + der_vy(i+1,j))+(4.0d0/(3.0d0))*(v_right-v_left))/dx
        

        uface(1) = 0.5d0*(u_left + u_right)
        uface(2) = 0.5d0*(v_left + v_right)


        u_left   = temperature(i+0,j)+(1.0d0/2.0d0)*derT(i+0,j,1)
        u_right  = temperature(i+1,j)-(1.0d0/2.0d0)*derT(i+1,j,1)
       

       derTa(1) = (0.5d0*(derT(i,j,1) + derT(i+1,j,1))+(4.0d0/(3.0d0))*(u_right-u_left))/dx
           
        
	    delV = derua(1)+derva(2) ! ux+vy+wz
	    
	    ! Txx, Txy
	    Txx = 2.d0*mu_lam*( derua(1) - delV/3.d0 )
	    Txy = mu_lam*( derua(2) + derva(1) )
	    
	    ! Tyx, Tyy
	    Tyx = Txy
	    Tyy = 2.d0*mu_lam*( derva(2) - delV/3.d0 )

	    
	    Theta(1) = uface(1)*Txx + uface(2)*Txy + kappa*derTa(1)
	    Theta(2) = uface(1)*Tyx + uface(2)*Tyy + kappa*derTa(2)
	    
	    vfx(1) = 0.d0
	    vfx(2) = Txx*sf(1) + Txy*sf(2)
	    vfx(3) = Tyx*sf(1) + Tyy*sf(2)
	    vfx(4) = Theta(1)*sf(1) + Theta(2)*sf(2)


	    fluxvX(i,j,1:4) = vfx(1:4)
        
     enddo
    enddo

    ! flux on J faces
    sf(1)=0; sf(2)=1
    
!     do k = 0, NZ
    do i = 1, NX
    
    	do j = 0, NY


	    	u_left  = u_vel(i,j+0)+(1.0d0/2.0d0)*der_uy(i,j+0)
	        u_right = u_vel(i,j+1)-(1.0d0/2.0d0)*der_uy(i,j+1)

	        v_left  = v_vel(i,j+0)+(1.0d0/2.0d0)*der_vy(i,j+0)
	        v_right = v_vel(i,j+1)-(1.0d0/2.0d0)*der_vy(i,j+1)

		   derua(1) = (0.5d0*(der_ux(i,j) + der_ux(i,j+1))+(4.0d0/(3.0d0))*(u_right-u_left))/dy
		   derva(1) = (0.5d0*(der_vx(i,j) + der_vx(i,j+1))+(4.0d0/(3.0d0))*(v_right-v_left))/dy
	     

	       derua(2) = (0.5d0*(der_uy(i,j) + der_uy(i,j+1))+(4.0d0/(3.0d0))*(u_right-u_left))/dy
	       derva(2) = (0.5d0*(der_vy(i,j) + der_vy(i,j+1))+(4.0d0/(3.0d0))*(v_right-v_left))/dy


	        uface(1) = 0.5d0*(u_left + u_right)
	        uface(2) = 0.5d0*(v_left + v_right)

	        u_left   = temperature(i,j+0)+(1.0d0/2.0d0)*derT(i,j+0,2)
	        u_right  = temperature(i,j+1)-(1.0d0/2.0d0)*derT(i,j+1,2)

	       	derTa(2) = (0.5d0*(derT(i,j,2) + derT(i,j+1,2))+(4.0d0/(3.0d0))*(u_right-u_left))/dy

        
		    delV = derua(1)+derva(2)
		    ! Txx, Txy
		    Txx = 2.d0*mu_lam*( derua(1) - delV/3.d0 )
		    Txy = mu_lam*( derua(2) + derva(1) )
		    
		    ! Tyx, Tyy
		    Tyx = Txy
		    Tyy = 2.d0*mu_lam*( derva(2) - delV/3.d0 )
		    
		    
		    Theta(1) = uface(1)*Txx + uface(2)*Txy + kappa*derTa(1)
		    Theta(2) = uface(1)*Tyx + uface(2)*Tyy + kappa*derTa(2)
		    
		    vfx(1) = 0.d0
		    vfx(2) = Txx*sf(1) + Txy*sf(2)
		    vfx(3) = Tyx*sf(1) + Tyy*sf(2)
		    vfx(4) = Theta(1)*sf(1) + Theta(2)*sf(2)


		    fluxvY(i,j,1:4) = vfx(1:4)
        
    	enddo
    enddo
       
        
    ! cells
    do j = 1, NY
   	 	do i = 1, NX
        
    		residual(i,j,1:4) = residual(i,j,1:4) + (fluxvX(i,j,1:4)-fluxvX(i-1,j,1:4))/dx + (fluxvY(i,j,1:4)-fluxvY(i,j-1,1:4))/dy
        
   		 enddo
    enddo
		



	end subroutine VF
	

	!***********************************************************************
	!*****                      Reconstruction		                   *****
	!***********************************************************************
	subroutine recon (NS, un, ulnew,urnew,mx,my,ghostp,n_eqn,reconstruction)

	implicit none
	integer				:: ix,NS, ghostp, n_eqn,i,j

	integer,intent(IN)  :: reconstruction
	double precision	:: un(-ghostp:NS+ghostp,n_eqn),prim(-ghostp:NS+ghostp,n_eqn)

	double precision	:: ulnew(-ghostp:NS+ghostp,n_eqn),urnew(-ghostp:NS+ghostp,n_eqn)
	double precision 	:: mx,my, lx, ly

	double precision,parameter :: kai=0.5d0

	lx = -my ; ly = mx;


	do ix =0,NS


	if(reconstruction .eq. 23)then
			!!!!!! 2nd and 3rd order reconstruction

			! 1 is density and central

			ulnew(ix,1) =        (1.0d0/2.0d0*un(ix,1)+1.0d0/2.0d0*un(ix+1,1))
			urnew(ix,1) =        ulnew(ix,1)

			! 4 is energy and central

			ulnew(ix,4) =        (1.0d0/2.0d0*un(ix,4)+1.0d0/2.0d0*un(ix+1,4))
			urnew(ix,4) =        ulnew(ix,4)

			
			! X-direction  rhov is central and rhou is upwind

			if(mx.eq.1) then

				ulnew(ix,3) =        (1.0d0/2.0d0*un(ix,3)+1.0d0/2.0d0*un(ix+1,3))
				urnew(ix,3) =        ulnew(ix,3)

				! ulnew(ix,2)   =un(ix  ,2)
			    ! urnew(ix,2)   =un(ix+1,2)

				ulnew(ix,2)   =-1.0d0/6.0d0*un(ix-1,2)+5.0d0/6.0d0*un(ix  ,2)+1.0d0/3.0d0*un(ix+1,2)
			    urnew(ix,2)   =+1.0d0/3.0d0*un(ix  ,2)+5.0d0/6.0d0*un(ix+1,2)-1.0d0/6.0d0*un(ix+2,2)

			endif

			! Y-direction  rhou is central and rhov is upwind

			if(my.eq.1) then

				ulnew(ix,2) =        (1.0d0/2.0d0*un(ix,2)+1.0d0/2.0d0*un(ix+1,2))
				urnew(ix,2) =        ulnew(ix,2)

				! ulnew(ix,3)   =un(ix  ,3)
			    ! urnew(ix,3)   =un(ix+1,3)
				ulnew(ix,3)   =-1.0d0/6.0d0*un(ix-1,3)+5.0d0/6.0d0*un(ix  ,3)+1.0d0/3.0d0*un(ix+1,3)
				urnew(ix,3)   = 1.0d0/3.0d0*un(ix  ,3)+5.0d0/6.0d0*un(ix+1,3)-1.0d0/6.0d0*un(ix+2,3)

			endif

	endif



	if(reconstruction .eq. 34)then
			!!!!!! 3rd and 4th order reconstruction

			! 1 is density and central

			ulnew(ix,1) =        kai*(-1.0d0/6.0d0*un(ix-1,1)+5.0d0/6.0d0*un(ix  ,1)+1.0d0/3.0d0*un(ix+1,1))+&
			             (1.0d0-kai)*(+1.0d0/3.0d0*un(ix  ,1)+5.0d0/6.0d0*un(ix+1,1)-1.0d0/6.0d0*un(ix+2,1))
			urnew(ix,1) =        ulnew(ix,1)

			! 4 is energy and central

			ulnew(ix,4) =        kai*(-1.0d0/6.0d0*un(ix-1,4)+5.0d0/6.0d0*un(ix  ,4)+1.0d0/3.0d0*un(ix+1,4))+&
			             (1.0d0-kai)*(+1.0d0/3.0d0*un(ix  ,4)+5.0d0/6.0d0*un(ix+1,4)-1.0d0/6.0d0*un(ix+2,4))
			urnew(ix,4) =        ulnew(ix,4)

			
			! X-direction  rhov is central and rhou is upwind

			if(mx.eq.1) then

				ulnew(ix,3) =        kai*(-1.0d0/6.0d0*un(ix-1,3)+5.0d0/6.0d0*un(ix  ,3)+1.0d0/3.0d0*un(ix+1,3))+&
				             (1.0d0-kai)*(+1.0d0/3.0d0*un(ix  ,3)+5.0d0/6.0d0*un(ix+1,3)-1.0d0/6.0d0*un(ix+2,3))
				urnew(ix,3) =        ulnew(ix,3)

				ulnew(ix,2)   =-1.0d0/6.0d0*un(ix-1,2)+5.0d0/6.0d0*un(ix  ,2)+1.0d0/3.0d0*un(ix+1,2)
			    urnew(ix,2)   =+1.0d0/3.0d0*un(ix  ,2)+5.0d0/6.0d0*un(ix+1,2)-1.0d0/6.0d0*un(ix+2,2)

			endif

			! Y-direction  rhou is central and rhov is upwind

			if(my.eq.1) then

				ulnew(ix,2) =        kai*(-1.0d0/6.0d0*un(ix-1,2)+5.0d0/6.0d0*un(ix  ,2)+1.0d0/3.0d0*un(ix+1,2))+&
				             (1.0d0-kai)*(+1.0d0/3.0d0*un(ix  ,2)+5.0d0/6.0d0*un(ix+1,2)-1.0d0/6.0d0*un(ix+2,2))
				urnew(ix,2) =        ulnew(ix,2)

				ulnew(ix,3)   =-1.0d0/6.0d0*un(ix-1,3)+5.0d0/6.0d0*un(ix  ,3)+1.0d0/3.0d0*un(ix+1,3)
				urnew(ix,3)   = 1.0d0/3.0d0*un(ix  ,3)+5.0d0/6.0d0*un(ix+1,3)-1.0d0/6.0d0*un(ix+2,3)

			endif

	endif

	if(reconstruction .eq. 56)then

			!!!!!! 5th and 6th order reconstruction

			! 1 is density and central

			ulnew(ix,1) =        kai*((1.0d0/60.0d0)*(+2.0d0*un(ix-2,1)-13.0d0*un(ix-1,1)+47.0d0*un(ix+0,1)+27.0d0*un(ix+1,1)-3.0d0*un(ix+2,1)))+&
			             (1.0d0-kai)*((1.0d0/60.0d0)*(-3.0d0*un(ix-1,1)+27.0d0*un(ix+0,1)+47.0d0*un(ix+1,1)-13.0d0*un(ix+2,1)+2.0d0*un(ix+3,1)))
			urnew(ix,1) =        ulnew(ix,1)


			! 4 is energy and central

			ulnew(ix,4) =        kai*((1.0d0/60.0d0)*(+2.0d0*un(ix-2,4)-13.0d0*un(ix-1,4)+47.0d0*un(ix+0,4)+27.0d0*un(ix+1,4)-3.0d0*un(ix+2,4)))+&
			             (1.0d0-kai)*((1.0d0/60.0d0)*(-3.0d0*un(ix-1,4)+27.0d0*un(ix+0,4)+47.0d0*un(ix+1,4)-13.0d0*un(ix+2,4)+2.0d0*un(ix+3,4)))
			urnew(ix,4) =        ulnew(ix,4)


			! X-direction  rhov is central (3) and rhou is upwind (2)
			if(mx.eq.1) then

				ulnew(ix,3) =        kai*((1.0d0/60.0d0)*(+2.0d0*un(ix-2,3)-13.0d0*un(ix-1,3)+47.0d0*un(ix+0,3)+27.0d0*un(ix+1,3)-3.0d0*un(ix+2,3)))+&
				             (1.0d0-kai)*((1.0d0/60.0d0)*(-3.0d0*un(ix-1,3)+27.0d0*un(ix+0,3)+47.0d0*un(ix+1,3)-13.0d0*un(ix+2,3)+2.0d0*un(ix+3,3)))
				urnew(ix,3) =        ulnew(ix,3)



			   ulnew(ix,2)   =(1.0d0/60.0d0)*(+2.0d0*un(ix-2,2)-13.0d0*un(ix-1,2)+47.0d0*un(ix+0,2)+27.0d0*un(ix+1,2)-3.0d0*un(ix+2,2))
		       urnew(ix,2)   =(1.0d0/60.0d0)*(-3.0d0*un(ix-1,2)+27.0d0*un(ix+0,2)+47.0d0*un(ix+1,2)-13.0d0*un(ix+2,2)+2.0d0*un(ix+3,2))



			endif

			! Y-direction  rhou is central (2) and rhov is upwind (3)

			if(my.eq.1) then

				ulnew(ix,2) =        kai*((1.0d0/60.0d0)*(+2.0d0*un(ix-2,2)-13.0d0*un(ix-1,2)+47.0d0*un(ix+0,2)+27.0d0*un(ix+1,2)-3.0d0*un(ix+2,2)))+&
				             (1.0d0-kai)*((1.0d0/60.0d0)*(-3.0d0*un(ix-1,2)+27.0d0*un(ix+0,2)+47.0d0*un(ix+1,2)-13.0d0*un(ix+2,2)+2.0d0*un(ix+3,2)))
				urnew(ix,2) =        ulnew(ix,2)


			    ulnew(ix,3)   =(1.0d0/60.0d0)*(+2.0d0*un(ix-2,3)-13.0d0*un(ix-1,3)+47.0d0*un(ix+0,3)+27.0d0*un(ix+1,3)-3.0d0*un(ix+2,3))
		        urnew(ix,3)   =(1.0d0/60.0d0)*(-3.0d0*un(ix-1,3)+27.0d0*un(ix+0,3)+47.0d0*un(ix+1,3)-13.0d0*un(ix+2,3)+2.0d0*un(ix+3,3))

			endif

	endif


	enddo
!*****************************!*****************************!*****************************!*****************************!*****************************

	
	end subroutine recon




 	