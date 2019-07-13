      program tov_one
      ! This program computes dimensionless pressure and energy density with a scaling factor E0 = 1.603E38
      ! actual pressure = E0*p
      ! gamma = 5/3
      ! EOS => P = k*(E)^(5/3)
      ! mass is computed in units of solar mass, m = actual mass/mo, where m0 is the solar mass.
      
      real :: r,rn,p,pn,h,k11,k21,k31,k41,t,k12,k22,k32,k42,m,mn,actp,actd,a
      real :: n
       write(*,*) "Enter the initial pressure in cgs units"
              read(*,*) a
              p = a/(1.613E38)
              n = 0
                

              open(1, file = "tov_sample.dat")
              r = 0.01
              h = 0.1
              m = 0.0
              
              do while(p > 0)

                            
              k11 = h*dp(r,p,m)
              k12 = h*dm(r,p)
              k21 = h*dp(r+h/2,p+k11/2,m + k12/2 )
              k22 = h*dm(r+h/2,p +k11/2 )
              k31 = h*dp(r+h/2,p+k21/2,m+k22/2)
              k32 = h*dm(r+h/2,p+k21/2)
              k41 = h*dp(r+h,p+k31,m+k32)
              k42 = h*dm(r+h,p+k31)

              pn = p +  (k11 + 2*k21 + 2*k31 + k41)/6                       
                            
              mn = m +  (k12 + 2*k22 + 2*k32 + k42)/6
              rn = r + h

              actp= p*1.613E38
              actd = (actp/5.38E9)**(3.0/5.0)

              write(1,*) actp,actd, m ,r

              p = pn
              r = rn
              m = mn
              n = n+1
              end do
              close(1)

              print*, "Radius of star in kilometres is:",r 
              !print*, "Mass of the star in units of solar mass is:", m
              !print*, n 

              

              

      
      end program tov_one

      function dp(r,p,m)
      implicit none
      real :: r,p,m
      real :: dp
      dp = -((p**(3.0/5.0)*(1 + 1.47*p**(2.0/5.0))*(m + 1.12*p*r**(3)))/((r**2) - 1.47*m*r))
      
      
      end function dp

      function dm(r,p)
      real p, r
      real dm
      dm = 0.7686*(p**(3.0/5.0))*(r**2)
      end function dm
