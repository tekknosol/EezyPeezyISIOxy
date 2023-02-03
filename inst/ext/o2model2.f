!======================================================================
! Bulk hypolimnetic oxygen model. 
! State variable: dissolved O2 concentration
! Units: m, g, day 
! Forcing variables: area, volume, temperature
!======================================================================

!======================================================================
! initializer for parameter common block
! N is number of parameters
!======================================================================
        subroutine siminit (odeparms)
            external odeparms
            integer,parameter :: N = 3
            double precision parms(N)
            common /myparms/parms
                call odeparms(N, parms)
            return
        end subroutine siminit

!======================================================================
! initializer for the forcing common block
! N is the number of forcing variables
!======================================================================
        subroutine simforc(odeforcs)
            external odeforcs
            integer,parameter :: N = 3
            double precision forcs(N)
            common /myforcs/forcs
                call odeforcs(N, forcs)
            return
        end subroutine simforc

!======================================================================
! Subroutine to calculate  derivatives 
! t: time, y: state, ydot: derivs, yout: output var
!======================================================================
        subroutine simderivs (neq, t, y, ydot, ip)                     ! yout, 

            implicit none

            double precision t, y, ydot                                ! , yout
            double precision Flux, Khalf, Theta          
            double precision Area_linear, Volume_linear, Temp_linear                    
            integer neq, ip(*)
            dimension y(neq), ydot(neq)                                 ! , yout(*)

            common /myparms /Flux, Khalf, Theta
            common /myforcs/Area_linear, Volume_linear, Temp_linear

!                if (ip(1) < 1) call rexit("nout should be at least 1")

! derivative ydot
                ydot = Area_linear * Flux *                             &
     &               y(1) / (y(1) + Khalf) *                            &
     &               Theta**(Temp_linear - 20.d0) /                     &
     &               Volume_linear

! output variables
!                yout(1) = ydot(1)
            return
        end subroutine simderivs
