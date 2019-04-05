      subroutine hypela2(d,g,e,de,s,t,dt,dum0,n,nn,kc,dum1,dum2,dum3,
     2                   disp,dispt,coord,ffn,frotn,strechn,eigvn,ffn1,
     3                   frotn1,strechn1,eigvn1,ncrd,itel,ndeg,ndm,
     4                   dum4,dum5,dum6,ifr,ifu)
c
c  
c     compute the 2-nd Piola-Kirchhoff tensor and the material tangent
c     for hyperelastic material
c     input    ffn1  current deformation gradient
c     output   s, d
c
c     creation date: 30.09.2016
c
      implicit real*8 (a-h,o-z)                                               dp
c
      include '../common/concom'
      include '../common/elmcom'
      include '../common/creeps'   ! common block in order to get the information about the time step
c
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      dimension xStr(3,3), xdTC(3,3,3,3)      
c
c     TIME STEP
c
c     xdt        \Delta t
c
c     MATERIAL PARAMETERS     
c
c     xk         \k
c     xmu        \mu
c     xKa        K
c     xm         m
c     xeta       \eta
c     xkappa     \kappa    
c     xbeta      \beta
c     xgamma     \gamma
c     xce        c
c      
c     INTERNAL VARIABLES
c     at time t={}^{n}t
c   
c     xnCi       C_i  
c     xnCii      C_ii
c     xns        s
c     xnsd       s_d
c
c     QUANTITIES WHICH REMAIN CONSTANT WITHIN THE TIMESTEP
c
c     xC         current Cauchy tensor          
c     xinvC      inverse of C
c     xdC        det(C)
c     integ      type of integration procedure
c                eq.1 MEBM (Modified Euler-Backwards Method)  
c                eq.2 EM   (Exponential Method)   
c     xP         deviatoric projection tensor (see eq. (28)) 
c     xU         second rank identity tensor
c
      dimension e(*),de(*),t(21),dt(21),g(*),d(ngens,ngens),s(ngens)
      dimension n(2),coord(ncrd,nnode),disp(ndeg,nnode),pstn1(3),
     2          dispt(ndeg,nnode),ffn(itel,itel),frotn(itel,itel),
     3          strechn(itel),eigvn(itel,itel),ffn1(itel,itel),
     4          frotn1(itel,itel),strechn1(itel),eigvn1(itel,itel),
     5          strain(3,3),einc(3,3),ffn1t(3,3),
     6          frotn1t(3,3),temp1(3,3), d1(6,6), s1(6)
c

      dimension xCi(3,3), xC1i(3,3), xC2i(3,3)
c      
c                 define material parameters
c
      data      xk/5.83333d0/,xmu/1.25d0/,xc1/0.000000001d0/,
     1          xc2/0.0000001d0/,xKa/0.00000001/,xm/1.0d0/,xeta/0.08d0/,   ! artificial viscosity for debugging!!!!!
     2          xkapa1/1.000d0/,xkapa2/1.00d0/,xbeta/1.00000d0/,
     3          xgamma/0.0000000001d0/
      interface
      
      
      
      function xprod (a,b)
      implicit real*8 (a-h,o-z)
      dimension xprod(3,3), a(3,3), b(3,3)
      end function xprod
      function xinver(a)
      implicit real*8 (a-h,o-z)
      dimension xinver(3,3), a(3,3)
      end function xinver
      function xsqrtr(a)
      implicit real*8 (a-h,o-z)
      dimension xsqrtr(3,3), a(3,3)
      end function xsqrtr
      function xten22(a,b)
      implicit real*8 (a-h,o-z)
      dimension xten22(3,3,3,3), a(3,3), b(3,3)
      end function xten22
      function xcros22(a,b)
      implicit real*8 (a-h,o-z)
      dimension xcros22(3,3,3,3), a(3,3), b(3,3)
      end function xcros22
      function xtransp(a)
      implicit real*8 (a-h,o-z)
      dimension xtransp(3,3), a(3,3)
      end function xtransp
      function xphi(a)
      implicit real*8 (a-h,o-z)
      dimension xphi(6), a(3,3)
      end function xphi
      function xinvph(a)
      implicit real*8 (a-h,o-z)
      dimension xinvph(3,3), a(6)
      end function xinvph
      function xpii(a)
      implicit real*8 (a-h,o-z)
      dimension xpii(6,6), a(3,3,3,3)
      end function xpii
      end interface   
c      open(777,file="parameters_set.txt",status='OLD')
c      read(UNIT=777,FMT=*) xgamma,xbeta,xc1,xc2,xkapa1,xkapa2
c      close(777)
c     write(UNIT=8888, FMT=*) xc1,xc2,xkapa1,xkapa2,xgamma,xbeta
c
c     computations for common block viskplast 
c
      
      xdt=timinc      
      xU=0.0d0                                        ! unit tensor 
      xU(1,1)=1.0d0
      xU(2,2)=1.0d0
      xU(3,3)=1.0d0
      xC=xprod(xtransp(ffn1),ffn1)                       ! current C
      xinvC=xinver(xC)                                   ! inverce of C
      xdC=xdet(xC)                                       ! det(C)
      xCiso=(xdC)**(-1.0d0/3.0d0)*xC                     ! isochoric part of C
      xnCiso=xprod(xtransp(ffn),ffn)                     ! previous C
      xnCiso=(xdet(xnCiso))**(-1.0d0/3.0d0)*xnCiso       ! isochoric part of {}^n C
c      
c     initial data for increment from array t
c     following storage convention is used:
c
c      t(1) not used (reserved for temperature)
c          xns=t(2)      
c          xnsd=t(3)
c          xnCi=xinvph(t(4:9))      
c          xnC1i=xinvph(t(10:15)) 
c          xnC2i=xinvph(t(16:21))  
c
c
      if (timinc.eq.0.0d0) then
c                          zero iteration      
                           t=0.0d0
			   t(4)=1.0d0
			   t(5)=1.0d0
			   t(6)=1.0d0 
			   t(10)=1.0d0
			   t(11)=1.0d0
			   t(12)=1.0d0
               t(16)=1.0d0
			   t(17)=1.0d0
			   t(18)=1.0d0                           
			   endif
c                   read current t	
      xns=t(2)
      xnsd=t(3)		   
      xnCi=xinvph(t(4:9))             ! Use MARC numeration to store a matrix
      xnC1i=xinvph(t(10:15))          ! Use MARC numeration to store a matrix
      xnC2i=xinvph(t(16:21))          ! Use MARC numeration to store a matrix
      xStSiz=dsqrt(xtrace(xprod(xCiso-xnCiso,xCiso-xnCiso)))  ! approx size of the strain increment
      if ((xdt.gt.0.0d0).and.(xStSiz.le.0.20d0))  then          
      call mainSub(xStr, xdTC, xCi, xC1i, xC2i, xs, xsd)
                                  else
				  xStr=0.0d0         ! zero iteration or the strain increment is too large
				  xdTC=0.0d0         ! zero iteration or the strain increment is too large
				  xCi=xnCi           ! no inelastic flow for zero iteration
				  xC1i=xnC1i         ! no inelastic flow for zero iteration
                  xC2i=xnC2i         ! no inelastic flow for zero iteration
				  xs=xns             ! no inelastic flow for zero iteration
				  xsd=xnsd           ! no inelastic flow for zero iteration
                                  endif
c
c     Update internal variables
c     Warning: only the increments of t should be defined
c     the array t is then updated automatically by MSC.MARC
c  
      dt(2)=xs-t(2)
      dt(3)=xsd-t(3)
      dt(4:9)=xphi(xCi)-t(4:9)
      dt(10:15)=xphi(xC1i)-t(10:15)
      dt(16:21)=xphi(xC2i)-t(16:21)

c
c     Output of stresses and tangent
c       
      if (ngens.eq.6) then                    
c                               GENERAL 3D:   ngens eq 6      
      s=xphi(xStr)                            ! MARC representation of stresses
      d=2.0d0*xpii(xdTC)                      ! MARC representation of tangent
                      else
c                               Axisymmetric but not plane strain! :  ngens eq 4
      s1=xphi(xStr)                           ! MARC representation of stresses
      d1=2.0d0*xpii(xdTC)                     ! MARC representation of tangent
      do 100 ii=1,4
      s(ii)=s1(ii)
      do 100 jj=1,4     
      d(ii,jj)=d1(ii,jj)
100   continue        

                      endif		      
      return
 1000 format (9e12.6)
      end ! end of hypela2 subroutine

c______________________________________________________________

      subroutine mainSub(xStr, xdTC, xCi1, xC1i1, xC2i1, xs1, xsd1)
c
c     subroutine to compute the second Piola-Kirchhoff stress
c                                    and the material tangent      
c     input:     C, nCi, nCi, ns, nsd   through the common block
c
c     output:        xStr      the stress tensor (2nd Piola-Kirchhoff)
c                    xdTC      the tangent
c                    xCi1      new value of Ci
c                   xCii1      new value of Cii 
c                     xs1      new value of s
c                    xsd1      new value of sd
c
      implicit real*8 (a-h,o-z)
      dimension xStr(3,3), xdTC(3,3,3,3), xCreal(3,3), 
     1          xStrPer(3,3), xPer(3,3), xCi(3,3), xC1i(3,3), xC2i(3,3),
     2          xCi1(3,3), xC1i1(3,3), xC2i1(3,3), xdif(6), xtem(6,6)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      interface
      function xStre(xCi)
      implicit real*8 (a-h,o-z)
      dimension xStre(3,3), xCi(3,3) 
      end function xStre
      function xnvph2(a)
      implicit real*8 (a-h,o-z)
      dimension xnvph2(3,3), a(6)
      end function xnvph2
      function xbas(i)
      implicit real*8 (a-h,o-z)
      dimension xbas(6)
      end function xbas      
      function xphi2(a)
      implicit real*8 (a-h,o-z)
      dimension xphi2(6), a(3,3)
      end function xphi2
      function xinver(a)
      implicit real*8 (a-h,o-z)
      dimension xinver(3,3), a(3,3)
      end function xinver
      function xinvpi(a)
      implicit real*8 (a-h,o-z)
      dimension xinvpi(3,3,3,3), a(6,6)      
      end function xinvpi
      end interface  
cccccccccccccccccccccccccccccccccccccccc         
      xCreal=xC                            ! save the real value of C for later use
c      
c         finite difference approximation
c
      xdelt=3.0d-8                                       ! increment used for numerical differentiation
      call SlHRobus(xksi, xCi1, xC1i1, xC2i1, xs1, xsd1,0,xksireuse) ! C is unpertubed     
      xStr=xStre(xCi1)                                   ! compute the unpertubed (real) stress tensor xStr
      xksireuse = xksi                                   ! save computed ksi for re-use 

      do 100 l=1,6
      xPer=xnvph2(xbas(l))*xdelt      
      xC=xCreal+xPer                                     ! Pertube C
      xinvC=xinver(xC)                                   ! update inverce of C
      xdC=xdet(xC)                                       ! update det(C)
      xCiso=(xdC)**(-1.0d0/3.0d0)*xC                     ! update isochoric part of C      
      call SlHRobus(xksi, xCi, xC1i, xC2i, xs, xsd,1,xksireuse)      ! C is pertubed, kreuse=2 means without update in D
      xStrPer=xStre(xCi)                                 ! compute the pertubed stress tensor
      xdif=xphi2(xStrPer-xStr)
      do 100 k=1,6
      xtem(k,l)=xdif(k)/xdelt
100   continue
      xdTC=xinvpi(xtem)
      end subroutine mainSub
      
c______________________________________________________________

      function xStre(xCi)
c
c     compute the second Piola-Kirchhoff stress
c     input: xCi
c     output xStre
c
      implicit real*8 (a-h,o-z)
      dimension xStre(3,3), xCi(3,3)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)

      interface
      function xprod (a,b)
      implicit real*8 (a-h,o-z)
      dimension xprod(3,3), a(3,3), b(3,3)
      end function xprod
      function xdev(a)
      implicit real*8 (a-h,o-z)
      dimension xdev(3,3), a(3,3)      
      end function xdev
      function xinver(a)
      implicit real*8 (a-h,o-z)
      dimension xinver(3,3), a(3,3)
      end function xinver 
      end interface
         ! the ansatz proposed by Hartmann and Neff (2003) is used here for the volumetric part  
      xJ=dsqrt(xdC)
      xStre=0.1d0*xk*(xJ**5.0d0-xJ**(-5.0d0) )*xinvC    ! volumetric part
      xStre=xStre+xmu*xprod(xinvC,xdev(xprod(xCiso,xinver(xCi))))
      end function xStre  
     
c______________________________________________________________

      subroutine SlHRobus(xksi, xCi, xC1i, xC2i, xs, xsd,kre,xksireuse)
c         Robust and efficient version of SolHD
c         based on the paper by Shutov (2015)
c         splitted Euler backward method is implemented
c
c         input:   C, nCi, nC1i, nC2i, ns, nsd --- through common block
c
c         output:  xksi, xCi, xC1i, xC2i, xs, xsd (current internal variables)
c
c         kre = 0 if unperturbed C isused
c         kre = 1 if C is perturbed, previously computed xksireuse can be re-used
c     
c
c
      implicit real*8 (a-h,o-z)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      dimension xCi(3,3), xC1i(3,3), xC2i(3,3), 
     1          xCiper(3,3), xC1iper(3,3), xC2iper(3,3)
     
c      interface
c      function xprod (a,b)
c      implicit real*8 (a-h,o-z)
c      dimension xprod(3,3), a(3,3), b(3,3)
c      end function xprod
c      end interface
c
c     compute trial driving force and trial isotropic hardening
c
      xFtr=xFF(xnCi, xnC1i, xnC2i)             ! trial driving force
      xtR=xgamma*(xns-xnsd)                    ! trial hardening
      if (xFtr.le.sqrt(2.0d0/3.0d0)*(xKa+xtR)) then
c
c                          the inelastic flow is frozen              
c
                       xksi=0.0d0
				       xCi=xnCi
				       xC1i=xnC1i
                       xC2i=xnC2i
				       xs=xns
				       xsd=xnsd
                                               return
                                               else
            ! use partionend Euler backward method
      xtoll=1.0d-13        ! convergence tolerance
      xdelt=1.0d-8         ! increment used for numerical differentiation
            if ((kre.eq.1).and.(xksireuse.ge.3.0d-8)) then
            xksi=xksireuse  ! re-use previously computed xksi
            call UpdateCiii(xCi, xC1i, xC2i, xksi, xtR)
            goto 80
            endif
       !  start with a single Newton-itartion for equation {H=0}
      xksi=0.0d0 ! zero approximation
      call UpdateCiii(xCiper, xC1iper, xC2iper, xksi+2.0d0*xdelt, xtR)
      xDvalP = xD(xksi+2.0d0*xdelt, xCiper, xC1iper, xC2iper)  ! perturbed value of H near ksi=0
      ! at this point xH(*,*,*,*) should not be used at all, since it may be ill-defined
      if (xDvalP.ge.0.0d0) then    ! xDvalP>0 means that the solution is between 0 and 2*xdelt
          ! find the solution by the bisection method
          kbisec=1                          ! kbisec eq 1 means that bisection is used
          xleft=0.0d0                       ! the left bound 
          xright=2.0d0*xdelt                ! the right bound 
          nbisec=0 ! number of bisection steps
70        continue
          nbisec=nbisec+1
          xcent=0.5d0*(xleft+xright)        ! bisection
       call UpdateCiii(xCiper, xC1iper, xC2iper, xcent, xtR)
          xDvalP = xD(xcent, xCiper, xC1iper, xC2iper)
          if (xDvalP.le.0.0d0) then
		      xleft=xcent
		                       else
              xright=xcent
		                endif
                        if (nbisec.ge.200) then
                        write(UNIT=8888, FMT=*)'Bisection went wrong' 
                        write(UNIT=8888, FMT=*)'xleft = ', xleft 
                        write(UNIT=8888, FMT=*)'xright = ', xright
                        pause 
                        endif
            if ((xright-xleft).gt.xtoll) goto 70
            xksi=0.5d0*(xleft+xright)              ! bisection finished
            nit=0
            call UpdateCiii(xCi, xC1i, xC2i, xksi, xtR)
            goto 300
                            endif
      kbisec=0  ! bisection not used
      xHval = xH(xksi, xnCi, xnC1i, xnC2i)  ! unperturbed value of H for ksi=0
      call UpdateCiii(xCiper, xC1iper, xC2iper, xksi+xdelt, xtR)
      xHvalP = xH(xksi + xdelt, xCiper, xC1iper, xC2iper)  ! perturbed value of H
      xdH = (xHvalP - xHval)/xdelt        ! derivative of H w.r.t. ksi
      xdksi = - xHval/xdH                 ! one Newton step
      xksi = xksi + xdksi
        if (xksi.lt.0.0d0) then
        write(UNIT=8888, FMT=*)'Negative ksi in {H=0} iteration' 
        write(UNIT=8888, FMT=*)'ksi=', xksi 
        end if
      call UpdateCiii(xCi, xC1i, xC2i, xksi, xtR)
      ! Newton iterations for {D=0}
80    nMax=15                             ! maximum number of iterations
      nit=0                               ! iteration counter
      xerr=xtoll+1.0d0                    ! current error
100   continue
      if ((xerr.le.xtoll).or.(nit.ge.nMax))  goto 300
      nit = nit +1
      xDval = xD(xksi, xCi, xC1i, xC2i)
      if (xerr.le.1.0d-8) goto 150 ! use previously computede derivative for small exerr
      call UpdateCiii(xCiper, xC1iper, xC2iper, xksi+xdelt, xtR)
      xDvalP = xD(xksi + xdelt, xCiper, xC1iper, xC2iper)  ! perturbed value of D
      xdD = (xDvalP - xDval)/xdelt        ! derivative of D w.r.t. ksi
           if (xdD.lt.0.0d0) then
           write(UNIT=8888, FMT=*)'Alaaaaaaaaaaaaaaaaaaaarm!'
           write(UNIT=8888, FMT=*)'Negative derivative of D!' 
           write(UNIT=8888, FMT=*)'Current value of C =', xC 
           write(UNIT=8888, FMT=*)'Current value of xksi =', xksi
           write(UNIT=8888, FMT=*)'kre =', kre
           write(UNIT=8888, FMT=*)'xksireuse =', xksireuse
           write(UNIT=8888, FMT=*)'xerr =', xerr
           write(UNIT=8888, FMT=*)'xDval =', xDval
           write(UNIT=8888, FMT=*)'xDvalP =', xDvalP
           endif
 150  continue
      xdksi = - xDval/xdD                 ! one Newton step
      xksi = xksi + xdksi
         if (xksi.lt.0.0d0) then
         write(UNIT=8888, FMT=*) 'xksi<0 in SlHDEf, 1st relax step',xksi  ! TEMP      
         write(UNIT=8888, FMT=*) 'nit =', nit  ! TEMP
         write(UNIT=8888, FMT=*) 'kbisec =', kbisec  ! TEMP  
         write(UNIT=8888, FMT=*) 'xHvalP, xDval =', xHvalP, xDval  ! TEMP     
         endif
      call UpdateCiii(xCi, xC1i, xC2i, xksi, xtR)
      xerr=abs(xdksi)
      goto 100
300   continue
          if (nit.ge.Nmax) then
          write(UNIT=8888, FMT=*) 'nit = Nmax! Error = ', xerr
          write(UNIT=8888, FMT=*) 'ksi = ', xksi
          write(UNIT=8888, FMT=*) 'kre = ', kre
          write(UNIT=8888, FMT=*) 'ksireuse = ', xksireuse
          write(UNIT=8888, FMT=*) 'current  C = ', xC
          write(UNIT=8888, FMT=*) 'xnCi = ', xnCi
          write(UNIT=8888, FMT=*) 'xnC1i = ', xnC1i
          write(UNIT=8888, FMT=*) 'xnC2i = ', xnC2i
          write(UNIT=8888, FMT=*) 'ns = ', ns
          write(UNIT=8888, FMT=*) 'nsd = ', nsd
          write(UNIT=8888, FMT=*) 'xHval = ', xHval
          write(UNIT=8888, FMT=*) 'xHvalP = ', xHvalP
          write(UNIT=8888, FMT=*) 'kbisec = ', kbisec
          write(UNIT=8888, FMT=*) 'xDval = ', xDval
          endif
       xs = xns + dsqrt(2.0d0/3.0d0)*xksi ! update s
       xsd = (xnsd + dsqrt(2.0d0/3.0d0)*xksi*xbeta*xs)/
     1       (1 + dsqrt(2.0d0/3.0d0)*xksi*xbeta)  
                                       return
                                       end if
       end subroutine SlHRobus 
         
c______________________________________________________________

      subroutine UpdateCiii(xCi, xC1i, xC2i, xksi, xtR)
c        
c         explicit update formulas + relaxation iterations
c
c         based on the paper by Shutov (2016) in CMAME
c       
c
c         input:   C     through common block
c                  xksi  as an argument
c                  xtR   as an argument
c
c         output:  xCi, xC1i, xC2i
c
c
      implicit real*8 (a-h,o-z)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      dimension xCi(3,3), xC1i(3,3), xC2i(3,3)
     
      interface
      function updtci(xksi, xesc1i, xesc2i, xtR, xfan)
      implicit real*8 (a-h,o-z)
      dimension updtci(3,3), xesc1i(3,3), xesc2i(3,3)
      end function updtci
      function xprod (a,b)
      implicit real*8 (a-h,o-z)
      dimension xprod(3,3), a(3,3), b(3,3)
      end function xprod
      end interface
c
c
c     initial approximation for Ci (as if without kinematic hardening)
c     0.005d0 is a scale factor used to reduce the stifnesses c1 and c2
c
      xCi = updtci(xksi, xnC1i, xnC2i, xtR, 0.005d0)  
      do 50 k=1,3  !relaxation interations pertaining to operator split
       xC1i =  xnC1i + xksi*xkapa1*xc1*xCi
       xC1i = (xdet(xC1i))**(-1.0d0/3.0d0)*xC1i
       xC2i =  xnC2i + xksi*xkapa2*xc2*xCi
       xC2i = (xdet(xC2i))**(-1.0d0/3.0d0)*xC2i
       xCi = updtci(xksi, xC1i, xC2i, xtR,1.0d0)
50    continue
       xC1i =  xnC1i + xksi*xkapa1*xc1*xCi
       xC1i = (xdet(xC1i))**(-1.0d0/3.0d0)*xC1i
       xC2i =  xnC2i + xksi*xkapa2*xc2*xCi
       xC2i = (xdet(xC2i))**(-1.0d0/3.0d0)*xC2i
       end subroutine UpdateCiii 

c______________________________________________________________
      
    
      function xD(xksi, xCi, xC1i, xC2i)
c
c        compute D=D(ksi) for ksi
c
c
c       input:  xksi       inelastic strain increment 
c      
      implicit real*8 (a-h,o-z)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      dimension xCi(3,3), xC1i(3,3), xC2i(3,3)            
      xtR=xgamma*(xns-xnsd)                               ! compute the trial hardening
      xR=xtR+sqrt(2.0d0/3.0d0)*xgamma*xksi                ! compute R(ksi)
      xR=xR/(1.0d0+sqrt(2.0d0/3.0d0)*xbeta*xksi)
      xD=(xeta*xksi/xdt)**(1.0d0/xm)
      xD=xD+dsqrt(2.0d0/3.0d0)*(xKa+xR)-xFF(xCi, xC1i, xC2i)
      end function xD
c______________________________________________________________


      function xH(xksi, xCi, xC1i, xC2i)
c
c        compute H=H(ksi) for ksi
c
c
c       input:  xksi       inelastic strain increment 
c      
      implicit real*8 (a-h,o-z)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      dimension xCi(3,3), xC1i(3,3), xC2i(3,3)            
      xtR=xgamma*(xns-xnsd)                               ! compute the trial hardening
      xR=xtR+sqrt(2.0d0/3.0d0)*xgamma*xksi                ! compute R(ksi)
      xR=xR/(1.0d0+sqrt(2.0d0/3.0d0)*xbeta*xksi)
      xH=xeta*xksi/xdt
      xH=xH-(xFF(xCi, xC1i, xC2i)-sqrt(2.0d0/3.0d0)*(xKa+xR))**xm
      end function xH
     
c     _____________________________________________________      


      function updtci(xksi, xesc1i, xesc2i, xtR,xfact)
      implicit real*8 (a-h,o-z)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      dimension updtci(3,3), xesc1i(3,3), xesc2i(3,3),
     1          xsqc1i(3,3), xsqc2i(3,3), xisq1i(3,3), xisq2i(3,3),
     2          xphi(3,3), xphi12(3,3), xAA(3,3), xY(3,3), xX(3,3)
      interface
      function xprod (a,b)
      implicit real*8 (a-h,o-z)
      dimension xprod(3,3), a(3,3), b(3,3)
      end function xprod
      function xinver(a)
      implicit real*8 (a-h,o-z)
      dimension xinver(3,3), a(3,3)
      end function xinver
      function xsqrtr(a)
      implicit real*8 (a-h,o-z)
      dimension xsqrtr(3,3), a(3,3)
      end function xsqrtr
      end interface
c
c     update the inelastic strain tensor Ci
c     input: xksi - ksi
c            xesc1i -  estimated {}^{n+1} C1i
c            xesc2i -  estimated {}^{n+1} C2i
c               xtR -  trial isotrpic hardening
c     all remaining varaibles - trough the common block
c     output updtci = updated Ci
c  
       xxc1=xc1*xfact   ! scale xc1 and xc2
       xxc2=xc2*xfact
       if (xksi.lt.0) then
       write(UNIT=8888, FMT=*)' Negative ksi in updtci', xksi
                      endif
       xR= (xtR + dsqrt(2.0d0/3.0d0)*xgamma*xksi)/
     1     (1.0d0+dsqrt(2.0d0/3.0d0)*xbeta*xksi)    ! current isotropic hardening
       xF = dsqrt(2.0d0/3.0d0)*(xKa+xR) + (xeta*xksi/xdt)**(1.0d0/xm) ! current overstress
c      xCiso is already available
       xcs= xxc1 + xxc2                                       ! not efficient, should be computed outside!!!
       xphi = (xxc1*xinver(xesc1i) + xxc2*xinver(xesc2i))/xcs ! Phi  not efficient, should be computed outside!!!
       xphi12 = xsqrtr(xphi)                                  ! Phi^(1/2)   not efficient, should be computed outside!!!
       xAA =  xnCi + (2.0d0*xksi*xmu/xF)*xCiso
       xAA = xprod(xphi12, xAA)
       xAA = xprod(xAA, xphi12)                       ! A is ready 
       xz0=(xdet(xAA)/xdet(xphi))**(1.0d0/3.0d0)      ! z0
       xz = xz0 - (xcs*xksi*xtrace(xAA))/(3.0d0*xz0*xF)  !
       xY = xz*xz*xU + 4.0d0*xcs*xksi/xF*xAA
       xY = xsqrtr(xY) + xz*xU
       xY = xinver(xY)
       xY = 2.0d0*xprod(xAA,xY)                       ! Y is ready
       xphi12 = xinver(xphi12) ! now xphi12 is Phi^(-1/2)
       xX = xprod(xphi12, xY)
       xX = xprod(xX, xphi12)
       updtci = (xdet(xX))**(-1.0d0/3.0d0)*xX
      end function updtci

c___________________________________________________________________ 

      function xFF(xCi, xC1i, xC2i)
c
c     compute the current driving force (combine eq. (114) and (121))
c
c     input      
c     xCi        C_i
c     xCii       C_ii
c
c     output
c     xFF        F           
c
      implicit real*8 (a-h,o-z)
      dimension xCi(3,3), xC1i(3,3), xC2i(3,3), xtemp(3,3)
      common /viscplast2/ xdt, xk, xmu, xKa, xm, xeta, xkapa1,
     1        xkapa2,xbeta, xgamma, xc1, xc2, xC, xdC, xCiso, xinvC, 
     2        xnCi, xnC1i, xnC2i, xns, xnsd, xU, xnCiso
      dimension xC(3,3), xCiso(3,3), xinvC(3,3), ! for the common block 
     1          xnCi(3,3), xnC1i(3,3), xnC2i(3,3), xU(3,3), xnCiso(3,3)
      interface
      function xprod (a,b)
      implicit real*8 (a-h,o-z)
      dimension xprod(3,3), a(3,3), b(3,3)
      end function xprod
      function xinver(a)
      implicit real*8 (a-h,o-z)
      dimension xinver(3,3), a(3,3)
      end function xinver
      function xdev(a)
      implicit real*8 (a-h,o-z)
      dimension xdev(3,3), a(3,3)      
      end function xdev
      end interface
      xtemp=xmu*xprod(xCiso,xinver(xCi))
      xtemp=xtemp-xc1/2.0d0*xprod(xCi,xinver(xC1i))
      xtemp=xtemp-xc2/2.0d0*xprod(xCi,xinver(xC2i))
      xtemp=xdev(xtemp)
      xFF=dsqrt(xtrace(xprod(xtemp,xtemp)))
      end function xFF

c___________________________________________________________________  

      function xprod(a,b)
      implicit real*8 (a-h,o-z)
      dimension xprod(3,3), a(3,3), b(3,3)
c
c     product of two 3x3 matrices
c     input: a, b
c     output xprod = a times b
c  
       do 100 i=1,3
       do 100 j=1,3
       xprod(i,j)=0.0d0
        do 100 k=1,3
        xprod(i,j)=xprod(i,j)+a(i,k)*b(k,j)
100   continue  
      end function xprod! end of xprod
      
c___________________________________________________________________  

      function xphi2(a)
      implicit real*8 (a-h,o-z)
      dimension xphi2(6), a(3,3), ii(6), jj(6)
c
c     convert second-rank tensor into a vector
c     input: a  (should be symmetric)
c     output xphi = OUR representation
c     
c     according to OUR convention 
c     3x3 symmetric matrix M is represented as a vector
c     (M11, M12, M13, M22, M23, M33) 
c      
      data ii/1,1,1,2,2,3/, jj/1,2,3,2,3,3/ !OUR numeration order
      do k=1,6
      xphi2(k)=a(ii(k),jj(k))
      end do     
      end function xphi2
      
c___________________________________________________________________        
      
      function xnvph2(a)
      implicit real*8 (a-h,o-z)
      dimension xnvph2(3,3), a(6), ii(6), jj(6)
c
c     inverse of xphi2
c      
      data ii/1,1,1,2,2,3/, jj/1,2,3,2,3,3/ !OUR numeration order
      do k=1,6
      xnvph2(ii(k),jj(k))=a(k)
      xnvph2(jj(k),ii(k))=a(k)
      end do     
      end function xnvph2
            
c___________________________________________________________________  
      function xphi(a)
      implicit real*8 (a-h,o-z)
      dimension xphi(6), a(3,3), ii(6), jj(6)
c
c     convert second-rank tensor into a vector
c     input: a  (should be symmetric)
c     output xphi = MARC representation
c     
c     according to MARC convention 
c     3x3 symmetric matrix M is represented as a vector
c     (M11, M22, M33, M12, M23, M31) 
c      
      data ii/1,2,3,1,2,3/, jj/1,2,3,2,3,1/ !MARC numeration order
      do k=1,6
      xphi(k)=a(ii(k),jj(k))
      end do     
      end function xphi
c___________________________________________________________________  

      function xbas(i)
      implicit real*8 (a-h,o-z)
      dimension xbas(6)
c
c     compute the standart basis in R6
c      
      xbas=0.0d0
      xbas(i)=1.0d0
      end function xbas
c___________________________________________________________________  

      function xinvph(a)
      implicit real*8 (a-h,o-z)
      dimension xinvph(3,3), a(6), ii(6), jj(6)
c
c     inverse of xphi
c      
      data ii/1,2,3,1,2,3/, jj/1,2,3,2,3,1/ !MARC numeration order
      do k=1,6
      xinvph(ii(k),jj(k))=a(k)
      xinvph(jj(k),ii(k))=a(k)
      end do     
      end function xinvph
 
c____________________________________________________________________ 

      function xinvpi(a)
      implicit real*8 (a-h,o-z)
      dimension xinvpi(3,3,3,3), a(6,6)
      dimension ii(6), jj(6)
c
c     convert a 6x6 matrix into a fourth-rank tensor
c     input: a
c     output xinvpi 
c     
c     according to OUR convention 
c       
      data ii/1,1,1,2,2,3/, jj/1,2,3,2,3,3/ !OUR numeration order            
      do 100 k=1,6
      do 100 l=1,6
      xinvpi(ii(k),jj(k),ii(l),jj(l))=a(k,l)
      if (ii(l).ne.jj(l)) then
      xinvpi(ii(k),jj(k),ii(l),jj(l))=0.5d0*a(k,l)          !!!!! IS IT TRUE ?????? 
                          endif
      xinvpi(jj(k),ii(k),ii(l),jj(l))=xinvpi(ii(k),jj(k),ii(l),jj(l))
      xinvpi(ii(k),jj(k),jj(l),ii(l))=xinvpi(ii(k),jj(k),ii(l),jj(l))
      xinvpi(jj(k),ii(k),jj(l),ii(l))=xinvpi(ii(k),jj(k),ii(l),jj(l))     
100   continue  
      end function xinvpi
 
c____________________________________________________________________  

      function xpii(a)
      implicit real*8 (a-h,o-z)
      dimension xpii(6,6), a(3,3,3,3)
      dimension ii(6), jj(6)
c
c     convert fourth-rank tensor into a 6x6 matrix
c     input: a
c     output xpii = representation material matrix
c     
c     according to MSC.MARC convention 
c     3x3 symmetric stress is represented as a vector
c     (S11, S22, S33, S12, S23, S31) 
c     
c     the strain tensor is probably represented by
c     (E11, E22, E33, 2*E12, 2*E23, 2*E31)
c       
      data ii/1,2,3,1,2,3/, jj/1,2,3,2,3,1/   !MARC numeration order            
      do 100 m=1,6
      do 100 n=1,6
       xpii(m,n)=0.5d0*a(ii(m),jj(m),ii(n),jj(n))
       xpii(m,n)=xpii(m,n)+0.5d0*a(ii(m),jj(m),jj(n),ii(n))      
c       if (ii(n).ne.jj(n)) xpii(m,n)=1.0d0*xpii(m,n)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100   continue  
      end function xpii
 
c____________________________________________________________________     
            
      function xten22(a,b)
      implicit real*8 (a-h,o-z)
      dimension xten22(3,3,3,3), a(3,3), b(3,3)
c
c     tensor product of two sdecond-rank tensors (in cartesian coordinates)
c     input: a, b
c     output xten22 = a \otimes b (see eq. (64)) 
c  
       do 100 i=1,3
       do 100 j=1,3
       do 100 k=1,3
       do 100 l=1,3
       xten22(i,j,k,l)=a(i,j)*b(k,l)
100   continue  
      end function xten22
 
c____________________________________________________________________ 
                 
      function xcros22(a,b)
      implicit real*8 (a-h,o-z)
      dimension xcros22(3,3,3,3), a(3,3), b(3,3)
c
c     cross-product of two sdecond-rank tensors (in cartesian coordinates)
c     input: a, b
c     output xcros22 = a \times b (see eq. (65)) 
c  
       do 100 i=1,3
       do 100 j=1,3
       do 100 k=1,3
       do 100 l=1,3
       xcros22(i,j,k,l)=a(i,k)*b(l,j)
100   continue  
      end function xcros22
 

c____________________________________________________________________  


      function xSc222(a,b)
      implicit real*8 (a-h,o-z)
      dimension xSc222(12), a(12,12), b(12)
c
c     product of a 12 x 12 matrix
c     with a vector
c     input: a, b
c     output xSc222 = a \cdot b
c  
       xSc222=0.0d0
       do 100 i=1,12      
       do 100 j=1,12
       xSc222(i)=xSc222(i)+a(i,j)*b(j)       
100   continue  
      end function xSc222
 
      
c____________________________________________________________________   
          
 
      function xsqrtr(a)
      implicit real*8 (a-h,o-z)
      dimension xsqrtr(3,3), a(3,3), xI(3,3), xB(3,3)
      data      xpi/3.141592653589793239d0/

      interface
      function xprod (a,b)
      implicit real*8 (a-h,o-z)
      dimension xprod(3,3), a(3,3), b(3,3)
      end function xprod
      end interface
c
c     square root of a symmetric positive definite 3x3 matrix
c     input: a
c     output xsqrttrnsr = a^(1/2)
c  
c    First, find the eigenvalues of a, using Smith's (1961) algorithm

      xp1 = a(1,2)**2.0d0 + a(1,3)**2.0d0 + a(2,3)**2.0d0
      if (xp1.eq.0.0d0) then 
c     a is diagonal.
       xeig1 = a(1,1)
       xeig2 = a(2,2)
       xeig3 = a(3,3)
                       else
      xq = (a(1,1)+a(2,2)+a(3,3))/3.0d0
      xp2=(a(1,1)-xq)**2.0d0 + (a(2,2)-xq)**2.0d0 + (a(3,3)-xq)**2.0d0
     1                                                      +2.0d0*xp1
      xp=dsqrt(xp2/6.0d0)
      xI=a
      xI(1,1)=xI(1,1)-xq
      xI(2,2)=xI(2,2)-xq
      xI(3,3)=xI(3,3)-xq
      xB = (1.0d0/xp)*xI      !xI is the 3x3 identity matrix
      xr = xdet(xB)/2.0d0
c         In exact arithmetic for a symmetric matrix  -1 <= r <= 1
c      but computation error can leave it slightly outside this range.
        if (xr.le.-1.0d0) then
        xphi = xpi/3.0d0
        else if (xr.ge.1.0d0)  then
        xphi = 0.0d0
        else
        xphi = acos(xr)/3.0d0
        endif
c          the eigenvalues satisfy eig3 <= eig2 <= eig1
      xeig1 = xq + 2.0d0*xp*cos(xphi)
      xeig3 = xq + 2.0d0*xp*cos(xphi + (2.0d0*xpi/3.0d0))
      xeig2 = 3.0d0*xq - xeig1 - xeig3 ! since trace(A) = eig1 + eig2 + eig3
                      endif
c       eigenvalues of the square root
      xeig1=xeig1**0.5d0
      xeig2=xeig2**0.5d0
      xeig3=xeig3**0.5d0
c       invariants of the square root
      xI1= xeig1 + xeig2 + xeig3
      xI3= xeig1*xeig2*xeig3
      xI2=0.50d0*(xI1**2.0d0 -xeig1**2.0d0 -xeig2**2.0d0 -xeig3**2.0d0)
      xB = xprod(a,a)    ! now xB stores a^2
      xB = xB + (xI2-xI1**2.0d0)*a   ! now xB stores a^2+[I2-I1^2]a
      xB(1,1) = xB(1,1) - xI1*xI3    
      xB(2,2) = xB(2,2) - xI1*xI3     
      xB(3,3) = xB(3,3) - xI1*xI3    ! now xB stores a^2+[I2-I1^2]a - I1*I3*1
      xsqrtr = (1.0d0/(xI3-xI1*xI2))*xB
      end function xsqrtr

c____________________________________________________________________ 
                     
      function xinver(a)
      implicit real*8 (a-h,o-z)
      dimension xinver(3,3), a(3,3)
c
c     inverse of a 3x3 matrix
c     input: a
c     output xinver = a^(-1)
c  
      xinver=a
      call inv3x3(xinver,xinver,det,1)  ! MARC Utility routine inv3x3 is used !!!!
      end function xinver

c____________________________________________________________________  
           
      function xtransp(a)
      implicit real*8 (a-h,o-z)
      dimension xtransp(3,3), a(3,3)
c
c     transposition of a 3x3 matrix
c     input: a
c     output xtransp = a^(T)
c  
       do 100 i=1,3
       do 100 j=1,3
       xtransp(i,j)=a(j,i)
100   continue       
      end function xtransp


c____________________________________________________________________

      function xdev(a)
      implicit real*8 (a-h,o-z)
      dimension xdev(3,3), a(3,3)
c
c     deviatoric part of a 3x3 matrix
c     input: a
c     output xdev = a^(D)
c  
      xtr=-xtrace(a)/3.0d0
      xdev=a
      do 100 i=1,3 
      xdev(i,i)=xdev(i,i)+xtr 
100   continue       
      end function xdev


c____________________________________________________________________        
            
      function xdet(a)
      implicit real*8 (a-h,o-z)
      dimension a(3,3)
c
c     determinant of a 3x3 matrix
c     input: a
c     output xdet = det(a)
c 
      xx=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)
      xy=-a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
      xdet=xx+xy
      end function xdet

c______________________________________________________________________      
            
      function xtrace(a)
      implicit real*8 (a-h,o-z)
      dimension a(3,3)
c
c     trace of a 3x3 matrix
c     input: a
c     output xtrace = tr(a)
c 
      xtrace=a(1,1)+a(2,2)+a(3,3)
      end function xtrace