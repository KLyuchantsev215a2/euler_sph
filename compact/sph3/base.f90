program base
!spatial integration method  SPH
   ! use ogpf

    
    integer N,i,count_hole,count_section,k1,k2!the number of particles that sample the environment
    integer step!counter for time steps
    integer sqn
    real*8 :: rho_0, T, l,CFL!density, total calculation time ,the size of the side of the square, Courant number
    real*8 :: S,m,h!body area, mass of a single particle , smoothing radius
    real*8 :: nu,mu,cs_0,E,k,eta,damping !material constants
    real*8 :: dh!indent for calculating the derived kernel through finite differences
    real*8 :: dt,time_calculated!time step, time during calculation
    
    real*8, allocatable :: x(:,:)
    real*8 :: Mtest(3,3)
    real*8 :: Mtest1(3,3)
  !  real*8, allocatable :: xplot(:,:,:)
    real*8, allocatable :: x_init(:,:)
    real*8, allocatable :: table(:,:)
    real*8, allocatable :: v(:,:)
   
    real*8, allocatable :: W(:,:)
    real*8, allocatable :: Wper1(:,:)
    real*8, allocatable :: Wper2(:,:)
    real*8, allocatable :: Wper3(:,:)!tmp
    real*8, allocatable :: Wper4(:,:)!tmp
    real*8, allocatable :: nabla_W(:,:,:)
    real*8, allocatable :: nabla_W_0(:,:,:)
    
    real*8, allocatable :: Couchy(:,:,:)
    real*8, allocatable :: PK1(:,:,:)
    real*8, allocatable :: F(:,:,:)
    real*8, allocatable :: Ci(:,:,:)
    real*8, allocatable :: Ci_new(:,:,:)
    real*8, allocatable :: vol(:)
    
    real*8, allocatable :: acc(:,:)
    real*8, allocatable :: x_0(:,:),x_n_1(:,:),x_n_2(:,:),x_n_1_2(:,:),x_n_3_2(:,:)
    real*8, allocatable :: v_0_0(:,:),v_n_1(:,:),v_n_2(:,:),v_n_1_2(:,:),v_n_3_2(:,:)
    
    integer, allocatable :: index_hole(:)
    integer, allocatable :: index_section(:)
    
     interface
        function Compute_W (xi,xj,h)
            real*8 :: xi(2)
            real*8 :: xj(2)
            real*8 :: h 
            real*8 :: Compute_W
        end function Compute_W
      
        function det (M)
            real*8 :: M(3,3)
            real*8 :: det
        end function det
               
        function trace (M)
            real*8 :: M(3,3)
            real*8 :: trace
        end function trace
        
        function dev (M)
            real*8 :: M(3,3)
            real*8 :: dev(3,3)
        end function dev
      
     end interface
    
    open (unit=1, file="input41.txt", status='old',    &
             access='sequential', form='formatted', action='read' )
    open (unit=2, file="output_x.txt", action='write')
    open (unit=3, file="output_C.txt", action='write')
    
       
    
    read (1, 1100) rho_0, T,nu, mu, l, dh,CFL,N 
    write (*, 1113) rho_0, T,nu, mu, l, dh,CFL,N
    
    sqn=21
    !S=(1.25*0.6-3.14*0.25*0.25/4)
    S=1
    m=rho_0*S/N
    
    k=2.0*mu*(1.0+nu)/(3.0*(1.0-2.0*nu))
    damping=0!0.003
    eta=1.0/25.0
    E=9.0*k*mu/(3.0*k+mu)

    cs_0=sqrt((E+4.0/3.0*mu)/rho_0)
    h=1.4*sqrt(m/rho_0)
    dt=CFL*h/(cs_0)
    
    allocate(vol(N))
    allocate(x(2,N))
   ! allocate(xplot(2,N,int(T/dt)))
    allocate(x_init(2,N))
    allocate(v(2,N))
    allocate(table(N,30))
    
    allocate(acc(2,N))
    allocate(x_0(2,N),x_n_1(2,N),x_n_2(2,N),x_n_1_2(2,N),x_n_3_2(2,N))
    allocate(v_0_0(2,N),v_n_1(2,N),v_n_2(2,N),v_n_1_2(2,N),v_n_3_2(2,N))
    allocate(W(N,N))
    
    allocate(Wper1(N,N))
    allocate(Wper2(N,N))
    allocate(Wper3(N,N))!tmp
    allocate(Wper4(N,N))!tmp
    
    !allocate(nabla_W(2,N,N))
    allocate(nabla_W_0(2,N,N))
    
    allocate(F(2,2,N))
    allocate(Ci(2,2,N))
    allocate(Ci_new(3,3,N))
    allocate(Couchy(2,2,N))
    allocate(PK1(2,2,N))
    !allocate(Cochi(2,2,N))
   
    vol=m/rho_0
        
    do i=1,N
        read (1, 1110) a,x(1,i),x(2,i)
    enddo
    
      do i=1,N
        read (1, 1110) a,v(1,i),v(2,i)
      enddo
      
    
  !  do i=1,N !razrez
  !      
   !     if ((x(1,i)<=0.7) * (abs(x(2,i))<0.001)) then
   !          x(2,i)=x(2,i)+0.001
  !      end if
    
    !enddo
    
    call Create_Table(x,h+dh,table,N)
    
   !  count_hole=0
   ! count_section=0
    
  !  do i=1,N
         
       ! if(      (sqrt((x(1,i)-0.25)**2+(x(2,i)-0.325)**2))<(0.25/2+0.001)          ) then
        !        count_hole=count_hole+1
       ! end if
  !       if(      ((x(1,i)<0.2)*(x(2,i)>0.4)*(x(2,i)<0.85))        ) then
  !              count_hole=count_hole+1
   !     end if
        
   !     if ( (x(1,i)>0.7)*(x(2,i)<=0.001))     then
  !              count_section=count_section+1
  !      end if
        
  !  enddo
    
    allocate(index_hole(count_hole))
    allocate(index_section(count_section))
     
  !  k1=1
  !  k2=1
 !   do i=1,N
        
   !     if(      ((x(1,i)<0.2)*(x(2,i)>0.4)*(x(2,i)<0.85))        ) then! if(      (sqrt((x(1,i)-0.25)**2+(x(2,i)-0.325)**2))<(0.25/2+0.001)          ) then
   !             index_hole(k1)=i
  !              k1=k1+1
    !    end if
        
        
    !    if ( (x(1,i)>0.7)*(x(2,i)<=0.001))     then
   !             index_section(k2)=i                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
   !             k2=k2+1
   !     end if
        
   ! enddo
   
    !call plot_init(x,N,count_hole,count_section,index_section,index_hole)

  ! v=0
   x_init=x
    
   
   
   call Compute_nabla_W(x,h,vol,N,W,Wper1,Wper2,Wper3,Wper4,nabla_W_0,dh,table)!tmp
   call Compute_F(vol,x,x_init,nabla_W_0,N,F,table)
   Ci=F
   call  OneStepMaxwell(F,mu,k,eta,dt,Ci,N,Couchy,Ci_new,PK1)
   Ci(1:2,1:2,1:N)=Ci_new(1:2,1:2,1:N)
    
    do step=1,int(T/dt)
        x_0=x
        v_0_0=v
        call Compute_Acceleration(N,h,dh,rho_0,mu,k,eta,damping,vol,F,Couchy,PK1,x_0,x_init,v,nabla_W_0,nabla_W,W,Wper1,Wper2,Wper3,Wper4,acc,count_hole,count_section,index_section,index_hole,Ci,Ci_new,table)
       ! x_n_1=x_0+dt*v_0_0
      !  v_n_1=v_0_0+dt*acc
      !  call Compute_Acceleration(N,h,dh,rho_0,mu,k,eta,damping,vol,F,Couchy,PK1,x_n_1,x_init,v_n_1,nabla_W_0,nabla_W,W,Wper1,Wper2,Wper3,Wper4,acc,count_hole,count_section,index_section,index_hole,Ci,Ci_new,table)
      !  x_n_2=x_n_1+dt*v_n_1     
      !  v_n_2=v_n_1+dt*acc
       ! x_n_1_2=3.0/4.0*x_0+1.0/4.0*x_n_2
       ! v_n_1_2=3.0/4.0*v_0_0+1.0/4.0*v_n_2
       ! call Compute_Acceleration(N,h,dh,rho_0,mu,k,eta,damping,vol,F,Couchy,PK1,x_n_1_2,x_init,v_n_1_2,nabla_W_0,nabla_W,W,Wper1,Wper2,Wper3,Wper4,acc,count_hole,count_section,index_section,index_hole,Ci,Ci_new,table)
      !  x_n_3_2=x_n_1_2+dt*v_n_1_2
      !  v_n_3_2=v_n_1_2+dt*acc
      !  x=1.0/3.0*x_0+2.0/3.0*x_n_3_2
       ! v=1.0/3.0*v_0_0+2.0/3.0*v_n_3_2
        v=v+dt*acc
        x=x+dt*v
        
        
        call Compute_F(vol,x,x_init,nabla_W_0,N,F,table) 
        call  OneStepMaxwell(F,mu,k,eta,dt,Ci,N,Couchy,Ci_new,PK1)
        Ci(1:2,1:2,1:N)=Ci_new(1:2,1:2,1:N)
        
        time_calculated=(real(step)*dt)
        
        !xplot(1:2,1:N,step)=x
    
        write (2,1111) x(1,1681)-x_init(1,1681),x(2,1681)-x_init(2,1681),time_calculated
        write (3,1112) Couchy(1,2,841),Couchy(1,1,841),Couchy(2,2,841),time_calculated
     
    enddo
    
  
    
    pause
    
  !  call  plot(xplot,N,int(T/dt))
    
    
    deallocate(vol)
    deallocate(x)
    deallocate(x_init)
    deallocate(v)
    
    deallocate(acc)
    deallocate(x_0,x_n_1,x_n_2,x_n_1_2,x_n_3_2)
    deallocate(v_0_0,v_n_1,v_n_2,v_n_1_2,v_n_3_2)
    deallocate(W)
    deallocate(nabla_W)
    deallocate(nabla_W_0)
    deallocate(table)
    
     1100 format (7f10.6,1i4)
    1113 format ("Density "1f10.6,/,"Time "1f10.6,/,"Poisson's ratio " 1f10.6,/,"Shear modulus " 1f10.6,/,"Side of a square " 1f10.6,/,"For finite difference " 1f10.6,/,"CFL " 1f10.6,/,"Particle count " 1i4)
    1110 format (1i12,1f25.0,1f20.0)
    1111 format (3f10.6)
1112     format (4f10.6)
         
    !1100 format (7f10.6,1i3)
   ! 1113 format ("Density "1f10.6,/,"Time "1f10.6,/,"Poisson's ratio " 1f10.6,/,"Shear modulus " 1f10.6,/,"Side of a square " 1f10.6,/,"For finite difference " 1f10.6,/,"CFL " 1f10.6,/,"Particle count " 1i3)
   ! 1110 format (1f22.0,1f23.0)
   ! 1111 format (3f10.6)
   ! 1112 format (4f10.6)
    
    end program base
    
    
    function Compute_W(xi,xj,h)
        real*8::xi(2)
        real*8::xj(2)
        real*8::h
        
        real*8::r(2)
        real*8::q
        real*8::C
        real*8::KER
        KER=0
        r=xi-xj
        q=sqrt(r(1)*r(1)+r(2)*r(2))/h
        C=1.0/(3.14159265358979323846*h*h)

        if((q>=0)*(q<=1)) then
               KER=C*(10.0 / 7.0)*(1.0-3.0/2.0*q*q*(1.0-q/2.0))
        end if
    
        if ((q > 1) * (q <=2)) then
            KER = C*(10.0 / 7.0)*(1.0/4.0)*(2.0 - q)*(2.0 - q)*(2.0 - q)
        end if
        
    
    Compute_W=KER
    end function Compute_W
    
    function det (M)
         real*8 :: M(3,3)
         real*8 ::det
         det=(M(1,1)*M(2,2)-M(1,2)*M(2,1))*M(3,3)
         
    end function det
        
    function trace (M)
            real*8 :: M(3,3)
            trace=M(1,1)+M(2,2)+M(3,3)
    end function trace
    
  
    
    
    
