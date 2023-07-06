      ! This code calculates the distribution of MSD, for experimental SPT
      ! This code has be written by S.S. Ashwin
      ! Please send suggestions/comments etc to ss.ashwin@gmail.com 
      ! The units anre in pixel (experimental) units 1 pixel= 65nm  

          module param
            ! Number of nucleosomes tracked in data set 
            integer,parameter::nmax=1342!2317 
            ! Number of frames
            integer,parameter::nfrm=100!449            i
            ! time in units of frames =Number of frames 
            integer,parameter::tmax=100!449 
            ! Definition of PI
            real(kind=8),parameter::pi=4.0d0*ATAN(1.0d0)
            ! Tolerance for convergence
            real(kind=8),parameter::tol=0.0000001d0
            ! Tolerance for convergence
            real(kind=8),parameter::tol2=0.0000001d0
            ! ddf: bin size for msd, ddmax: bin size distance 
            real(kind=8)::ddf,ddmax
            ! maximum size of time mesh <= nframe
            integer,parameter::ntmax=100
            ! time at which we want to calculate P(M,tref)
            integer::tref
            ! keep trefmax > tref
            integer,parameter::trefmax=30
            ! maximum size of the distance mesh
            integer,parameter::ndmax=10
            !number of bins of the MSD mesh
            integer,parameter::dmesh=500
            ! number of bins of the distance mesh  
            integer,parameter::irmax=500
         end module
           program pm1d 
               use param  
               ! x,y coordinates ; dx and dy differences 
               ! msd=mean square disp ; nrm7 =normalization
               real(kind=8),dimension(nmax,0:tmax)::x,y,dx,dy,msd,nrm7
               !maxtime: maximum time a nucleosome with index is tracked
               !dd1: dummy var used ; 
               real(kind=8),dimension(nmax)::maxtime,dd1
               ! intime: frame in which the nucl. appears first
               ! outtime: frame from which the nucl. exits observation
               integer,dimension(nmax)::intime,outtime
               ! variables defined when used 
               real(kind=8)::xo,yo,dxdx,dydy,dxdy,dydx,d2,nsum,tt,ddx,ddy
               real(kind=8)::dd0,nconv,fc,nsum1,nsum3,nsum2,nsum0,v2,rv1,rv2
               integer::i,j,ji,j0,i0,k,j1,j2,onep,k1,k2,kmax,k0,delt
               real(kind=8)::a,b,a1,a2,kr,a0,b0,r2,nf,ns
               real(kind=8)::dd,x0,y0
               integer::k3,k4,ial,isl,ifs,k5,iot,icheck,ki
               ! arrays used in MSD calculation
               real(kind=8),dimension(0:ntmax)::dr2,nrm2
               ! arrays used incalculation of  van Hove Correlation
               real(kind=8),dimension(0:irmax,0:ntmax)::Gs,Gsn,Gsno
               real(kind=8),dimension(0:irmax,0:ntmax)::Gss,Gsf
               real(kind=8),dimension(0:ntmax)::nrm
               ! Array used to store Gaussian function
               real(kind=8),dimension(0:irmax,dmesh,0:ntmax)::gxd
               ! Array used for MSD mesh in RL agorithm
               real(kind=8),dimension(dmesh)::Df
               ! Array used for P(M) in RL algorithm
               real(kind=8),dimension(dmesh,0:ntmax)::pd,pdn,pdo
               !Dmax: max of MSD mesh, Dmin :min MSD mesh
               real(kind=8)::Dmax,Dmin
               real(kind=8)::rsm,rfm,ro,rt
               ! kc bin index dividing fast and slow nucleosomes
               integer,dimension(0:10)::kc
               real(kind=8)::x2 
               ! velocity array and xy components
               real(kind=8),dimension(1:nmax,0:tmax)::vel,velx,vely
               ! Array consisting of max number of nucleosomes
               integer,dimension(0:10)::nmx 

              ! Input file
               OPEN(unit=21,file="./20170202c059_RPE1_CAG_H2B_Halotag_TMR80pM_nonlam_non_starve_ctrl_info.txt")
               print*,"provide tref < ",ntmax
               read(*,*)tref
                kc=0
                kc(1)=7
                nmx=0
                nmx(1)=nmax

          ! PARAMETERS for trial P(M) used in RL algo 
               a1=1.15!0.50d0
               a2=0.0!4.380d0
               a=1.0d0
               b=0.0!1.0d0
               a0=0.59!0.1d0
               b0=0.0!0.7d0    
          ! Define the MSD mesh for RL algo
                Dmax= 20.0d0!
                Dmin=0.0d0  
                ddf=(Dmax-Dmin)/real(dmesh)
          ! Initial ize arrays used for van Hove corrln
                Gs=0.0d0 
                Gsn=0.0d0 
                Gsno=0.0d0 
           ! Defining real space mesh bin
                ddmax=real(ndmax)/real(irmax) 
          
                nsum=0.0d0
                vel=0.0d0
                velx=0.0d0
                vely=0.0d0

        do ic=1,1
                ! This part of the code reads the data file
                ! containing SPT nucleosome trajectories
                ! It assigns xy coord to each nucleosome
                ! It assigns intime/outime the first frame/last frame
                ! each of the nucleosomes are seen in
                x=0
                y=0
                dd=0.0d0   
                i=0
                ji=0
                i0=0
                j0=1
                nrm=1.0d0 
                intime=0
                maxtime=0.0d0
                outtime=0
                do 
                    read(20+ic,*)j0,i,xo,yo
                     print*,j0,i,"nmax",nmax 
                     x(j0,i)=xo
                     y(j0,i)=yo
                     if(j0/=ji)then
                         intime(j0)=i
                         maxtime(ji)=real(i0)
                         x0=xo
                         y0=yo
                         if(ji/=0)then
                           dd1(ji)=dd0
                           if(dd0>tol)then
                           endif
                         endif
                         ji=j0
                         i0=1
                         dd0=0.0d0
                    else
                         i0=i0+1
                         if(i0>tmax)then
                             print*,"i0",i0,j0
                             stop
                        endif
                            outtime(j0)=i
                    endif
                      
                    if(j0==nmx(ic))then
                      maxtime(j0)=i0
                      exit
                    endif
                enddo
                ! Here we write out the nucleosome trajectory for each nucleosome
                do k=1,nmax
                      do j=intime(k),outtime(k)
                         write(60,*)x(k,j),y(k,j)
                      enddo   
                      write(60,*)
                enddo 
                ! Some nucleosome appear only once
                ! we assign intime=outtime for such nucleosomes
                ! calculate velocity of each nucle.
                do j=1,nmx(ic)
                   if(outtime(j)==0)then
                       outtime(j)=intime(j)
                   else
                       do k=intime(j)+1,outtime(j)
                          ddx=x(j,k)-x(j,k-1)
                          ddy=y(j,k)-y(j,k-1)
                           r2=ddx*ddx+ddy*ddy
                           r2=SQRT(r2)
                           vel(j,k-1)=r2
                           velx(j,k-1)=ddx
                           vely(j,k-1)=ddy
                        enddo
                   endif 
               enddo

               nrm2=1.0d0
               dr2=0.0d0
               msd=0.0d0
               nrm7=1.0d01
               ! Calculation of MSD 
               do j=1,nmx(ic)
                    do k1=intime(j),outtime(j)-1
                         do k2=k1,outtime(j)
                            ddx=x(j,k1)-x(j,k2)
                            ddy=y(j,k1)-y(j,k2)
                            r2=ddx*ddx+ddy*ddy
                            dr2(k2-k1) = dr2(k2-k1)+r2
                            nrm2(k2-k1) = nrm2(k2-k1)+1.0d0
                            msd(j,k2-k1)=msd(j,k2-k1) + r2
                            nrm7(j,k2-k1)=nrm7(j,k2-k1)+1.0d0
                         enddo
                    enddo 
               enddo
              
               do k1=1,ntmax
                     if(nrm2(k1)>tol)then
                        write(9,*)k1,dr2(k1)/nrm2(k1)!msd(j,k1)/nrm7(j,k1)
                     endif
               enddo
               print*,"msd done"
              
              ! Calculate experimental van Hove Correlation

               ifs=0
               isl=0
               Gss=0.0d0
               Gsf=0.0d0
               Gs=0.0d0
               nsum0=0.0d0
               nsum1=0.0d0
               nsum3=0.0d0 
               do j=1,nmx(ic)
                        delt=outtime(j)-intime(j)
                        do k1=intime(j),outtime(j)-1
                              do k2=k1+1,outtime(j)
                                     dxdx=x(j,k2)-x(j,k1)   
                                     dydy=y(j,k2)-y(j,k1)  
                                     r2=dxdx*dxdx+dydy*dydy
                                     r2=SQRT(r2) 
                                     k4=r2/(ddmax)
                                     k3=k2-k1 
                                     if(k4 <= irmax)then
                                           Gs(k4,k3)=Gs(k4,k3)+1.0d0
                                     endif
                            enddo
                       enddo
               enddo    
                   
           
               do i=1,ntmax
                    nsum=0.0d0 
                    do j=1,irmax
                         Gs(j,i)=Gs(j,i)/real(2.0d0*pi*j*ddmax) 
                         nsum=nsum + Gs(j,i)
                    enddo
                    nsum1=0.0d0
                    do j=1,irmax
                       if(nsum>tol)then
                          Gs(j,i)=Gs(j,i)/nsum
                          nsum1=nsum1+Gs(j,i)
                       endif
                    enddo
               enddo

               do i=1,irmax
                 write(363,*)real(i)*ddmax,Gs(i,tref)
               enddo 
            

 !  Start Richardson-Lucy 

              !  calc g(x|D,t)::g(x,D,t)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do k=1,dmesh
                   Df(k)=real(k)*ddf
                   do i=1,ntmax
                        do j=1,irmax
                           x2=real(j*j*ddmax*ddmax)
                           gxd(j,k,i)=(1.0d0/(pi*Df(k)))*exp(-x2/(Df(k)))  
                        enddo
                   enddo 
               enddo
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !  Trial P(M)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              pd=0.0d0
              nsum2=0.0d0
              do i=1,ntmax
                 do k=1,dmesh
                    kr=real(k)*ddf
                    pd(k,i)=a1*exp(-a*(kr-a0)**2.0d0)!+a2*exp(-b*(kr-b0)**2.0d0)
                 enddo
              enddo

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! Richardson-Lucy Iteration
             ! Gn and pd iteration
             !!!!!!!!!!!!!!!!!!:wq!!!!!!!!!!!!!!!!
             pdn=pd 
             Gsn=Gs
             do
                   pdo=pdn 
                   call integrate1(gxd,pdn,Gsn)
                   do i=1,trefmax
                      nsum=0.0d0
                      do j=1,irmax
                         nsum=nsum+Gsn(j,i)
                      enddo
                      do j=1,irmax
                         Gsn(j,i)=Gsn(j,i)/real(nsum)
                      enddo 
                   enddo
                   
                   Gsno=Gsn
                   call integrate2(gxd,Gs,Gsn,pdn) 
                   do i=1,trefmax
                      nsum1=0.0d0
                      do k=1,dmesh
                         nsum1=nsum1+pdn(k,i)
                      enddo
                      if(nsum1>tol)then 
                        do k=1,dmesh
                         pdn(k,i)=pdn(k,i)/nsum1
                       enddo
                      endif
                  enddo   
                  nconv=0.0d0
                  do k=1,dmesh
                      do i=1,trefmax
                         nconv=nconv+(pdn(k,i)-pdo(k,i))**2.0d0
                      enddo
                  enddo   
                  print*,"nconv",nconv,"tol2",tol2  
                  if(nconv<tol2)then
                     exit
                  endif

                  do k=1,dmesh
                      write(970,*)real(k*ddf), pdn(k,20)
                   enddo
                      write(970,*)  
                     flush(970) 
              enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
             do i=tref,tref
                nsum1=0.0d0 
                do k=1,dmesh
                   nsum1=nsum1+pdn(k,i)
                enddo
                do k=1,dmesh
                   pdn(k,i)=pdn(k,i)/real(nsum1)
                enddo
             enddo
             do i=tref,tref
                nsum1=0.0d0
                do k=1,dmesh
                   if((k==1).or.(k==dmesh))then
                     fc=0.50d0
                   else
                     fc=1.0d0
                   endif
                   nsum1=nsum1+pdn(k,i)*fc*ddf
                enddo
                do k=1,dmesh
                   pdn(k,i)=pdn(k,i)/nsum1
                enddo 
                nsum1=0.0d0
                do k=1,dmesh
                   if((k==1).or.(k==dmesh))then
                      fc=0.50d0
                   else
                      fc=1.0d0
                   endif
                   nsum1=nsum1+pdn(k,i)*fc*ddf
                enddo
               print*,"NORM",nsum1
            enddo  
            do i=1,ntmax!tref,tref
               do j=1,irmax
                  nsum1=0.0d0
                  do k=1,kc(ic)
                     nsum1=pdn(k,i)*gxd(j,k,i)*ddf + nsum1
                  enddo
                  Gss(j,i)=nsum1 
                  nsum1=0.0d0
                  do k=kc(ic)+1,dmesh
                     nsum1=pdn(k,i)*gxd(j,k,i)*ddf + nsum1
                  enddo
                  Gsf(j,i)=nsum1 
               enddo
             enddo
             ns=0.0d0
             nf=0.0d0
             do i=tref,tref
                do k=1,kc(ic)
                   ns=ns+pdn(k,i)
                   write(2001,*)k,pdn(k,i)
                enddo
                do k=kc(ic)+1,dmesh
                   nf=nf+pdn(k,i)
                   write(2002,*)k,pdn(k,i)
                enddo
             enddo  
             do i=1,trefmax
                  do k=1,dmesh
                     write(750+ic,*)real(k)*ddf,pdn(k,i)!*ifs/(isl+ifs)
                     write(770+ic,*)real(k)*ddf*65.0d0*65.0d0,pdn(k,i)!*ifs/(isl+ifs)
                  enddo
                  write(750+ic,*)
                  write(770+ic,*)
             enddo       
             do i=tref,tref
                  do j=1,irmax
                       peri=real(2.0d0*pi*j*real(ddmax)*65.0)
                       write(851,*)real(j*real(ddmax)*65.0),Gsf(j,i)*peri
                       write(852,*)real(j*real(ddmax)*65.0),Gss(j,i)*peri
                       write(853,*)real(j*real(ddmax)*65.0),Gsn(j,i)*peri,Gs(j,i)*peri
                  enddo
                  write(851,*)
                  write(852,*)
                  write(853,*)
            enddo   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !         enddo      !ic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo  !ic      
   end program   

             


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
      subroutine  integrate2(gxd,Gs,Gsn,pdn)
              use param
              real(kind=8),dimension(1:irmax,dmesh,1:ntmax),intent(in)::gxd
              real(kind=8),dimension(1:irmax,1:ntmax),intent(in)::Gs
              real(kind=8),dimension(dmesh,1:ntmax),intent(inout)::pdn
              real(kind=8),dimension(1:irmax,1:ntmax),intent(in)::Gsn
              real(kind=8)::nsum1,nsum2,fc
              integer::i,j,k    
              
              ! do i=1,ntmax
               do i=1,trefmax
                  do k=1,dmesh
                    nsum1=0.0d0
                    do j=1,irmax
                        if((j==1).or.(j==irmax))then
                          fc=0.50d0
                        else
                          fc=1.0d0
                        endif  
                       if(Gsn(j,i)>tol)then
                          nsum1=nsum1 + fc*Gs(j,i)*gxd(j,k,i)*real(ddmax)*(real(2.0d0*pi*j*ddmax))/Gsn(j,i) 
                       else
                          nsum1=nsum1 + fc*gxd(j,k,i)*real(ddmax)*(real(2.0d0*pi*j*ddmax)) 
                       endif
                    enddo
                    pdn(k,i)=nsum1*pdn(k,i)
                enddo
              enddo
       end subroutine      
                   
                  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine  integrate1(gxd,pd,Gsn)
              use param
              real(kind=8),dimension(1:irmax,dmesh,1:ntmax),intent(in)::gxd
              real(kind=8),dimension(dmesh,1:ntmax),intent(in)::pd
              real(kind=8),dimension(1:irmax,1:ntmax),intent(out)::Gsn
              real(kind=8)::nsum1,nsum2,fc
              integer::i,j,k   

               Gsn=0.0d0          
                   ! do i=1,ntmax
                    do i=1,trefmax
                         do j=1,irmax
                             nsum1=0.0d0
                             do k=1,dmesh
                                if((k==0).or.(k==dmesh))then
                                   fc=0.50d0
                                else
                                   fc=1.0d0
                                endif  
                                nsum1 = nsum1 + fc*gxd(j,k,i)*pd(k,i)*ddf
                             enddo 
                             Gsn(j,i)=nsum1
                        enddo
                    enddo
       end subroutine    
   
