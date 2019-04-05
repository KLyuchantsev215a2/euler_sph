subroutine Compute_nabla_W(x,h,vol,N,W,Wper1,Wper2,Wper3,Wper4,nabla_W,dh,table)!tmp
        integer N
        real*8::x(2,N)
        real*8::xper1(2,N)
        real*8::xper2(2,N)
        real*8::xper3(2,N)!tmp
        real*8::xper4(2,N)!tmp
        real*8::h
        real*8::retur
        real*8::vol(N)
        real*8::W(N,N)
        real*8::Wper1(N,N)
        real*8::Wper2(N,N)
        real*8::Wper3(N,N)!tmp
        real*8::Wper4(N,N)!tmp
        real*8::nabla_W(2,N,N)
        real*8:: dh
        real*8 :: table(N,30)
        xper1=x
        xper2=x
        xper3=x!tmp
        xper4=x!tmp
        xper1(1,1:N)=xper1(1,1:N)+dh
        xper2(2,1:N)=xper2(2,1:N)+dh
        xper3(1,1:N)=xper3(1,1:N)-dh!tmp
        xper4(2,1:N)=xper4(2,1:N)-dh!tmp
        
        !call Compute_W_cor(x,x,h,N,vol,W)
        call Compute_W_cor(x,xper1,h,N,vol,Wper1,table)
        call Compute_W_cor(x,xper2,h,N,vol,Wper2,table)
        call Compute_W_cor(x,xper3,h,N,vol,Wper3,table)!tmp
        call Compute_W_cor(x,xper4,h,N,vol,Wper4,table)!tmp
        
         nabla_W(1,1:N,1:N)=(Wper1(1:N,1:N)-Wper3(1:N,1:N))/(2*dh)
         nabla_W(2,1:N,1:N)=(Wper2(1:N,1:N)-Wper4(1:N,1:N))/(2*dh)
        !nabla_W(1,1:N,1:N)=(Wper1(1:N,1:N)-W(1:N,1:N))/dh 
        !nabla_W(2,1:N,1:N)=(Wper2(1:N,1:N)-W(1:N,1:N))/dh  
        
        return
    end 