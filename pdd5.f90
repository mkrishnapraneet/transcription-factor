          module param
            ! Number of nucleosomes tracked in data set 
            integer,parameter::nmax=1096!2317 
            ! Number of frames
            integer,parameter::nfrm=250!449            i
           ! time in units of frames =Number of frames 
            integer,parameter::tmax=250!449 
          
            integer,parameter::ilmax=10000
            ! Definition of PI
            real(kind=8),parameter::pi=4.0d0*ATAN(1.0d0)
            ! Tolerance for convergence
            real(kind=8),parameter::tol=0.0000001d0
            ! Tolerance for convergence
            real(kind=8),parameter::tol2=0.0000001d0

            real(kind=8),parameter::alpha=0.5d0
            real(kind=8)::ddf,ddmax
            integer,parameter::ntmax=100
            integer,parameter::ivcmax=100
            integer::tref
            integer,parameter::trefmax=30
            integer,parameter::tref0=1
            integer,parameter::dref=10
            integer,parameter::ndmax=10
            integer,parameter::ndf=500
            integer,parameter::dmesh=500
            real(kind=8),parameter::rll=2.0d0   
            integer,parameter::irmax=500
            real(kind=8),parameter::convt=0.04545
            real(kind=8),parameter::convd=0.065
            real(kind=8),parameter::wmax= 5.0d0!2.0d0*pi
            real(kind=8)::dw
            real(kind=8),parameter::ccrit=0.740d0
         end module
           program matrixdata
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
               real(kind=8)::xo,yo,dxdx,dydy,dxdy,dydx,d2,nsum,tt,ddx,ddy
               real(kind=8),parameter::deltd=0.05d0
               integer::i,j,ji,j0,i0,k,j1,j2,onep,k1,k2,kmax,k0,delt
               real(kind=8)::sumdiag,sumnondiag,dt,fact,dd,x0,y0
               real(kind=8),dimension(ilmax)::distd  
               !real(kind=8),dimension(1:irmax,0:tref-1)::vel,velx,vely
             !  real(kind=8),dimension(0:tref)::dos
               real(kind=8),dimension(0:ntmax)::dos,dr2,nrm2,vhist,dr2s,dr2f,nrm1,nrm0,dr2s0
               real(kind=8),dimension(0:ntmax)::drsl,drfs,nrmsl,nrmfs,vxsl,vxfs,vysl,vyfs
               real(kind=8)::vx2,vy2,r200
               real(kind=8),dimension(0:irmax)::velcor,velcorsl,velcorfs,vass,nrm4,vaff,nrm5,vch
               real(kind=8),dimension(nmax,0:irmax)::vassx,vassy
             !  real(kind=8),dimension(0:tref-1)::velcor,velcorsl,velcorfs
               real(kind=8),dimension(0:irmax,0:ntmax)::Gs,Gsn,Gsno,G2dto1d
               real(kind=8),dimension(0:irmax,0:ntmax)::Gfc,Gsc
               real(kind=8),dimension(irmax)::GR
               real(kind=8)::FRQ,FIQ,Q,FRQ1
               real(kind=8),dimension(0:irmax,0:ntmax)::Gsf,Gss
               real(kind=8),dimension(0:irmax)::Gt
               real(kind=8),dimension(0:ntmax)::nrm
               real(kind=8),dimension(0:irmax,dmesh,0:ntmax)::gxd
               real(kind=8),dimension(dmesh)::Df
               real(kind=8),dimension(dmesh,0:ntmax)::pd,pdn,pdo
               real(kind=8)::Dmax,Dmin,dd0,nconv,fc,nsum1,nsum3,nsum2,nsum0,v2,rv1,rv2
               integer,dimension(nmax)::cat,trajcat
               real(kind=8)::rsm,rfm,ro,rt,msdcrit
               integer,dimension(0:10)::kc
               real(kind=8)::phi,t0,dtmax,x2,relsf,relss,relff 
               real(kind=8),dimension(1:nmax,0:tmax)::vel,velx,vely
               real(kind=8),dimension(0:irmax)::vrsl,vrfs,vcss,vcsall,avcss,cvcss,gcss,gcfs,gcsall,gcfall,gcall
               real(kind=8),dimension(0:irmax)::vcfall,vcall
               real(kind=8),dimension(0:irmax)::avcff,cvcff,vcff
               real(kind=8),dimension(0:irmax)::avcfs,cvcfs,gcff,vcfs
               real(kind=8),dimension(0:irmax)::avcsall,cvcsall
               real(kind=8),dimension(0:irmax)::avcfall,cvcfall
               real(kind=8),dimension(0:irmax)::avcall,cvcall
               real(kind=8),dimension(0:irmax)::avgGss,avgGff,avgGfs,avgGsall,avgGfall,avgGall
               real(kind=8),dimension(0:irmax)::avgvcss,avgvcff,avgvcfs,avgvcsall,avgvcfall,avgvcall
               real(kind=8),dimension(0:irmax)::avgcvcss,avgcvcff,avgcvcfs,avgcvcsall,avgcvcfall,avgcvcall
             !  real(kind=8),dimension(0:tref-1)::velcor,velcorsl,velcorfs
               !real(kind=8),dimension(1:ntmax)::nsum2 
               real(kind=8)::a,b,a1,a2,kr,a0,b0,r2,avgr2,avgdf,nf,ns,thst,dv,vccc,avgr
               integer::k3,k4,ial,isl,ifs,k5,itraj,itot,iot,icheck,ki
               real(kind=8)::gammass,gammasf,gammaff,delta1,delta2,delta3,delta4,peri
               integer,dimension(0:10)::nmx 
               real(kind=8),dimension(0:ivcmax,0:irmax)::pvf,pvs
               real(kind=8),dimension(nmax,0:irmax)::nrm6
               real(kind=8),dimension(1:irmax,1:irmax)::G2D,G2Dn,G2Dno
               real(kind=8),dimension(dmesh,dmesh)::P2D,P2Dn,P2Do,cnorm
               real(kind=8),dimension(1:irmax,dmesh)::gxd2D
               real(kind=8),dimension(1:irmax)::G2Dj1norm,G2Dj2norm
               OPEN(unit=22,file="41_RPE1_CAG_H2B.txt")
       
               read(*,*)tref
               kc(0)=29!24
               kc(1)=37
               kc(2)=7
               kc(3)=119
               kc(4)=30
               kc(5)=28
               kc(6)=17
               kc(7)=15
               kc(8)=35
               kc(9)=20
               kc(10)=30

               nmx(0)=3614
               nmx(1)=3024
               nmx(2)=nmax
               nmx(3)=2185
               nmx(4)=2618
               nmx(5)=3739
               nmx(6)=2801
               nmx(7)=3750
               nmx(8)=3795
               nmx(9)=2732
               nmx(10)=3265


          ! PARAMETERS 
               a1=1.15!0.50d0
               a2=0.0!4.380d0
               a=1.0d0
               b=0.0!1.0d0
               a0=0.59!0.1d0
               b0=0.0!0.7d0    
               dos=0.0d0                
                dw=1/real(tref)!wmax/real(tref)
                dv=1.0d0/real(ivcmax)
                Dmax= 20.0d0!0.50d0
                Dmin=0.0d0  
                ddf=(Dmax-Dmin)/real(dmesh)
                Gs=0.0d0 
                Gsn=0.0d0 
                Gsno=0.0d0 
                ddmax=real(ndmax)/real(irmax) 
                nsum=0.0d0
                itraj=0
                vel=0.0d0
                velx=0.0d0
                vely=0.0d0
                itot=0
                pvs=0.0d0
                pvf=0.0d0

               do ic=2,2
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
                  !  if(j0==nfrm)then
                      maxtime(j0)=i0
                      exit
                    endif
              enddo
       !   enddo  
          do k=1,nmax
           do j=intime(k),outtime(k)
              write(60,*)x(k,j),y(k,j)
           enddo   
           write(60,*)
          enddo 
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
             nrm7=1.0d0
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
                   if(k1==tref)then
                      msdcrit=dr2(k1)/nrm2(k1)
                   endif     
                 endif
           enddo


            ifs=0
            isl=0
            Gss=0.0d0
            Gsf=0.0d0
            Gs=0.0d0
            nsum0=0.0d0
            nsum1=0.0d0
            nsum3=0.0d0 
            Gt=0.0d0
            k5=int(2.0d0/ddmax)
            cat=0 
            print*,"starting calc" 

              do j=1,nmx(ic)
                        delt=outtime(j)-intime(j)
                        do k1=intime(j),outtime(j)-1
                              do k2=k1+1,outtime(j)
                                     dxdx=x(j,k2)-x(j,k1)  !x(j,k2)-x(j,k1)
                                     dydy=y(j,k2)-y(j,k1)  !x(j,k2)-x(j,k1)
                                     r2=dxdx*dxdx+dydy*dydy
                                     r2=SQRT(r2) 
                                     k4=r2/(ddmax)
                                     k3=k2-k1!ABS(k2-k1) 
                                     if(k4 <= irmax)then
                                                Gs(k4,k3)=Gs(k4,k3)+1.0d0
                                    endif
                            enddo
                       enddo
               enddo    
                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


           
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
            

!
              !  calc g(x|D,t)::g(x,D,t)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do k=1,dmesh
                Df(k)=real(k)*ddf
             do i=1,ntmax
                do j=1,irmax
                   x2=real(j*j*ddmax*ddmax)
                   !gxd(j,k,i)=1.0d0/(4.0d0*pi*Df(k)*(real(i)**alpha))*exp(-x2/(4.0d0*pi*Df(k)*(real(i)**alpha)))  
                   gxd(j,k,i)=(1.0d0/(4.0d0*pi*Df(k)))*exp(-x2/(Df(k)))  
                enddo
             enddo 
           enddo
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  Trial P(D)
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
         ! Gn and pd iteration
         !!!!!!!!!!!!!!!!!!:wq!!!!!!!!!!!!!!!!
          pdn=pd 
         ! Gsn=Gs
         Gsn=Gs
          do
                   pdo=pdn 
                   call integrate1(gxd,pdn,Gsn)
                  ! do i=1,ntmax
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
                   !Gs=Gsn
                   call integrate2(gxd,Gs,Gsn,pdn) 
                  ! do i=1,ntmax
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
                     ! do i=1,ntmax
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
           !do i=1,ntmax!tref,tref!ntmax
           do i=tref,tref!ntmax
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
          stop
          print*,"starting G2D"
!   G2D 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   G2D trial <\delta(|r(tref0)-r(0)|-r1)\delta(|r(tref)-r(0)|-r2)>
             G2D=0.0d0   
              do j=1,nmx(ic)
                        delt=outtime(j)-intime(j)
                      if(delt>=tref)then
                        do k1=intime(j),outtime(j)-tref
                                     dxdx=x(j,k1+tref)-x(j,k1)  !x(j,k2)-x(j,k1)
                                     dydy=y(j,k1+tref)-y(j,k1)  !x(j,k2)-x(j,k1)
                                     r2=dxdx*dxdx+dydy*dydy
                                     r2=SQRT(r2) 
                                     k4=r2/(ddmax)
                                     print*,"r2",r2,x(j,k1+tref),x(j,k1)
                                     print*,"k3",k3
                                     dxdx=x(j,k1+tref0)-x(j,k1)  !x(j,k2)-x(j,k1)
                                     dydy=y(j,k1+tref0)-y(j,k1)  !x(j,k2)-x(j,k1)
                                     r2=dxdx*dxdx+dydy*dydy
                                     r2=SQRT(r2) 
                                     k3=r2/(ddmax) 
                                     if((k4< irmax).and.(k3<irmax))then
                                                G2D(k3,k4)=G2D(k3,k4)+1.0d0
                                    endif
                        enddo
                  endif
              enddo   

              do j1=1,irmax
                 G2Dj1norm=0.0
                 do j2=1,irmax
                    G2Dj1norm(j1)=G2Dj1norm(j1) + G2D(j1,j2)
                 enddo
             enddo
              do j2=1,irmax
                 G2Dj2norm=0.0
                 do j1=1,irmax
                    G2Dj2norm(j2)=G2Dj2norm(j2) + G2D(j1,j2)
                 enddo
             enddo
             do j1=1,irmax
               do j2=1,irmax
                 G2D(j1,j2)=G2D(j1,j2)/(4.0d0*pi*j1*j2*ddmax*ddmax)
               enddo
             enddo  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gxd2D ! 2D Gaussian distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
           do k1=1,dmesh
              !  k1=100
                Df(k1)=real(k1)*ddf
                do j1=1,irmax
                   x2=real(j1*j1*ddmax*ddmax)
                   gxd2D(j1,k1)=exp(-x2/(Df(k1)))!(1.0d0/(4.0d0*pi*Df(k1)))*exp(-x2/(Df(k1))) 
                enddo
           enddo
           cnorm=1.0d0
           do k1=1,dmesh
              do k2=1,dmesh
                 do j1=1,irmax
                     do j2=1,irmax
                       cnorm(k1,k2)=cnorm(k1,k2)+4.0d0*pi*pi*ddmax*j1*ddmax*j2*gxd2D(j1,k1)*gxd2D(j2,k2)*ddmax*ddmax
                      enddo
                 enddo
                  !  cnorm(k1,k2)=cnorm(k1,k2)
              enddo  
           enddo     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! P2D trial
          print*,"starting p2D trial"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
         do k1=1,dmesh
            do k2=1,dmesh
               P2D(k1,k2)=pdn(k1,tref0)*pdn(k2,tref)
            enddo
         enddo   
        do j1=1,dmesh
           do j2=1,dmesh
              if(P2D(j1,j2)>tol)then
               write(6300,*)real(ddf*j1),real(ddf*j2),-Log(P2D(j1,j2))
              endif
           enddo
        enddo
         flush(6300) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   nsum=0.0d0
                   do j2=1,irmax
                      do j1=1,irmax
                         nsum=nsum+G2D(j1,j2)
                      enddo
                   enddo   
                      do j1=1,irmax
                        do j2=1,irmax
                         G2D(j1,j2)=G2D(j1,j2)/real(nsum)
                        enddo 
                      enddo 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !iteration
          print*,"starting iteration"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! G2D and P2D iteration
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         P2Dn=P2D 
         G2Dn=G2D
          do
                   P2Do=P2Dn 
                   print*,"start integ12D"
                   call integrate12D(gxd2D,cnorm,P2Dn,G2Dn)
                   print*,"end integ12D"
                      nsum=0.0d0
                   do j2=1,irmax
                      do j1=1,irmax
                         nsum=nsum+G2Dn(j1,j2)
                      enddo
                   enddo   
                      do j2=1,irmax
                        nsum1=0.0d0
                        nsum2=0.0d0
                        do j1=1,irmax
                         G2Dn(j1,j2)=G2Dn(j1,j2)/real(nsum)
                         nsum1=nsum1 + G2Dn(j1,j2)
                         nsum2=nsum2 + G2D(j1,j2)
                        enddo 
                   !     write(8000+j,*)real(j2*ddmax),nsum1,nsum2
                      enddo 
                   !   flush(8000+j)
                   
                   print*,"nsum",nsum
                   G2Dno=G2Dn
                   !Gs=Gsn
                   print*,"start integ22D"
                   call integrate22D(gxd2D,cnorm,G2D,G2Dn,P2Dn) 
                   print*,"end integ22D"
                   nsum1=0.0d0
                   do k2=1,dmesh
                      nsum=0.0d0
                      do k1=1,dmesh
                         nsum1=nsum1+P2Dn(k1,k2)
                         nsum=nsum + P2Dn(k1,k2)
                      enddo
                    !     write(7000+j,*)real(k2*ddf),nsum
                   enddo   
                  ! flush(7000+j)
                   do k2=1,dmesh
                      do k1=1,dmesh
                         P2Dn(k1,k2)=P2Dn(k1,k2)/real(nsum1)
                      enddo
                   enddo   
                   print*,"nsum1",nsum1
                   nconv=0.0d0
                   do k1=1,dmesh
                      do k2=1,dmesh
                         nconv=nconv+(P2Dn(k1,k2)-P2Do(k1,k2))**2.0d0
                      enddo
                   enddo   
                   
                   print*,"nconv",nconv,"tol2",tol2  
                   if(nconv<tol2)then
                     exit
                   endif

           enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print*,"end iteration"
         do k1=1,dmesh
            do k2=1,dmesh
             if(P2Dn(k1,k2)>tol)then
             write(6500,*) real(k1*ddf),real(k2*ddf),-Log(P2Dn(k1,k2))
             endif
            enddo
         enddo    
         do k2=1,dmesh
            nsum1=0.0d0
            do k1=1,dmesh
             !if(P2Dn(k1,k2)>tol)then
                nsum1=nsum1+P2Dn(k1,k2)      
            ! endif
            enddo
            write(6504,*)real(k2*ddf),nsum1
         enddo 

                    do j1=1,irmax
                       nsum1=0.0d0
                        do j2=1,irmax
                         nsum1=G2Dn(j1,j2)+ nsum1
                        enddo 
                         peri=real(2.0d0*pi*j1*real(ddmax))
                        write(6501,*)real(j1*ddmax),nsum1*peri,Gsn(j1,tref0)*peri
                    enddo 

                    do j2=1,irmax
                       nsum1=0.0d0
                       nsum=0.0d0
                        do j1=1,irmax
                         nsum1=G2Dn(j1,j2)+ nsum1
                         nsum=G2D(j1,j2) + nsum
                        enddo 
                         peri=real(2.0d0*pi*j2*real(ddmax))
                        write(6502,*)real(j2*ddmax),nsum1*peri,Gsn(j2,tref)*peri,nsum*peri
                    enddo 

                   do j1=1,irmax
                        do j2=1,irmax
                        if((G2Dn(j1,j2)>tol).and.(G2D(j1,j2)>tol))then 
                         peri=4.0d0*pi*pi*real(j1*j2)*ddmax*ddmax       
                         write(6600,*)real(j1*ddmax),real(j2*ddmax),(G2Dn(j1,j2)*peri)
                         write(6400,*)real(j1*ddmax),real(j2*ddmax),(G2D(j1,j2)*peri)
                        endif
                    enddo 
                   enddo  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      enddo  !ic      
         end program   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             


          subroutine  integrate22D(gxd2D,cnorm,G2D,G2Dn,P2Dn)
              use param
              real(kind=8),dimension(1:irmax,dmesh),intent(in)::gxd2D
              real(kind=8),dimension(dmesh,dmesh),intent(in)::cnorm
              real(kind=8),dimension(1:irmax,1:irmax),intent(in)::G2D
              real(kind=8),dimension(dmesh,dmesh),intent(inout)::P2Dn
              real(kind=8),dimension(1:irmax,1:irmax),intent(in)::G2Dn
              real(kind=8)::nsum1,nsum2,fc1,fc2,term1
              integer::i,j,k1,k2,j1,j2    
              
               do k1=1,dmesh
                  do k2=1,dmesh
                    nsum1=0.0d0
                    do j1=1,irmax
                       do j2=1,irmax
                        if((j1==1).or.(j1==irmax))then
                          fc1=0.50d0
                        else
                          fc1=1.0d0
                        endif  
                        if((j2==1).or.(j2==irmax))then
                          fc2=0.50d0
                        else
                          fc2=1.0d0
                        endif  
                      term1=gxd2D(j1,k1)*gxd2D(j2,k2)*(real(4.0d0*pi*j1*j2*ddmax*ddmax))*ddmax/cnorm(k1,k2) 
                    !   term1=1.0d0
                       if(G2Dn(j1,j2)>tol)then
                          nsum1=nsum1+ term1*(G2D(j1,j2)/G2Dn(j1,j2)) 
                       else
                          nsum1=nsum1 + term1 
                       endif
                     enddo
                    enddo 
                    P2Dn(k1,k2)=nsum1*P2Dn(k1,k2)
                enddo
              enddo
       end subroutine      
                   
                  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine  integrate12D(gxd2D,cnorm,P2D,G2Dn)
              use param
              real(kind=8),dimension(1:irmax,dmesh),intent(in)::gxd2D
               real(kind=8),dimension(dmesh,dmesh),intent(in)::cnorm
              real(kind=8),dimension(dmesh,dmesh),intent(in)::P2D
              real(kind=8),dimension(1:irmax,1:irmax),intent(out)::G2Dn
              real(kind=8)::nsum1,nsum2,fc1,fc2
              integer::i,j1,k1,j2,k2   

               G2Dn=0.0d0          
                    do j1=1,irmax
                         do j2=1,irmax
                             nsum1=0.0d0
                             do k1=1,dmesh 
                              do k2=1,dmesh
                                if((k1==0).or.(k1==dmesh))then
                                   fc1=0.50d0
                                else
                                   fc1=1.0d0
                                endif  
                                if((k2==0).or.(k2==dmesh))then
                                   fc2=0.50d0
                                else
                                   fc2=1.0d0
                                endif  
                                nsum1 = nsum1 + fc1*fc2*gxd2D(j1,k1)*P2D(k1,k2)*gxd2D(j2,k2)*ddf*ddf/cnorm(k1,k2)
                              enddo 
                             enddo 
                             G2Dn(j1,j2)=nsum1
                        enddo
                    enddo

       end subroutine    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      subroutine fourier(Q,FR,GR,GI)
                 use param
                 INTEGER :: I,J,INM
                 REAL(KIND=8),INTENT(IN) ::Q 
                 REAL(KIND=8), INTENT (IN), DIMENSION (irmax) :: FR
                 REAL(KIND=8), INTENT (OUT) :: GR,GI
                 REAL(KIND=8)::R,dtheta,theta,nsum1,nsum2
                 REAL(KIND=8)::fc1,fc2

                 INM=500
                 dtheta=2.0d0*pi/real(INM)
                 nsum1=0.0d0
                 nsum2=0.0d0
                 DO I=1,irmax
                    R=ddmax*real(I)
                    if((I==0).or.(I==irmax))then
                      fc2=1.0d0
                    else
                      fc2=0.5d0
                    endif
                    DO J=1,INM
                       if((J==0).and.(J==INM))then
                         fc1=1.0d0
                       else
                         fc1=0.50d0
                       endif
                       theta=real(J)*dtheta
                       nsum1=nsum1 + FR(I)*fc1*fc2*COS(-Q*R*COS(theta))*ddmax*R*dtheta 
                       nsum2=nsum2 + FR(I)*fc1*fc2*SIN(-Q*R*COS(theta))*ddmax*R*dtheta 
                     !  print*,COS(Q*R*COS(theta)),theta,COS(theta)
                    ENDDO
               ENDDO
               GR=nsum1/pi
               GI=msum2/pi
         end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          subroutine phonondos(vc,dos)
                use param   
               ! real(kind=8),dimension(0:tref),intent(in)::vc
                real(kind=8),dimension(0:ntmax),intent(in)::vc
              !  real(kind=8),dimension(0:tref),intent(out)::dos
                real(kind=8),dimension(0:ntmax),intent(out)::dos
                real(kind=8)::zw
                integer::i,j
                 dos=0.0d0
                 do i=0,ntmax
                    w=real(i)*dw
                    zw=0.0d0
                    do j=0,tref-1
                       zw=zw+2.0d0*COS(w*real(j))*vc(j)/vc(0) 
                    enddo
                    dos(i)=zw
                enddo
           end subroutine  
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
   
