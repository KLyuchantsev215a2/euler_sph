subroutine OneStepMaxwell(F,mu,k,eta,dt,Ci,N,Couchy,Ci_new,PK1)
    
    integer:: N
    real*8 :: F(2,2,N)
    real*8 :: Couchy(2,2,N)
    real*8 :: Ci(2,2,N)
    real*8 :: Ci_new(3,3,N)
    real*8 :: mu
    real*8 :: dt
    real*8 :: k
    real*8 :: eta
    
    real*8 ::PK1(2,2,N)
    real*8 ::PK1_trans_tmp(2,2)
    
    real*8::detFp
    real*8 ::Fp(3,3)
    real*8 ::invFp(3,3)
    real*8 :: trans_Fp(3,3)
    
    real*8 ::C(3,3)
    real*8 ::invC(3,3)
    real*8 ::Ci3x3(3,3)
    
    real*8 ::Cip(3,3)
    real*8 ::invCip(3,3)
    real*8 ::multCCi(3,3)
    
    real*8 ::C_iso(3,3)
    real*8 ::dev_C_iso(3,3)
    real*8 ::Stress2PK(3,3)
    real*8::devmultCCi(3,3)
    real*8:: Couchy_tmp(3,3)
    real*8:: Couchy_tmp1(3,3)
    
    Couchy=0
    PK1=0
    
    do i=1,N
        Ci3x3=0
        C=0
        Fp=0
        
  
        Fp(1:2,1:2)=F(1:2,1:2,i)
        
        Fp(3,3)=1
        
        detFp=(Fp(1,1)*Fp(2,2)-Fp(1,2)*Fp(2,1))
        
        call trans(Fp,trans_Fp)
          
        call mymulty(trans_Fp,Fp,C)

        C_iso=detFp**(-2.0/3.0)*C   ! isochoric part of current C

        Ci3x3(1:2,1:2)=Ci(1:2,1:2,i)
        Ci3x3(3,3)=1.0/(Ci(1,1,i)*Ci(2,2,i)-Ci(1,2,i)*Ci(2,1,i))
        
        Cip=Ci3x3;   ! C_i from the previous time step at the current point
        
        Ci3x3 = Ci3x3 + dt*mu/eta*C_iso
        
        detCi3x3=(Ci3x3(1,1)*Ci3x3(2,2)-Ci3x3(1,2)*Ci3x3(2,1))*Ci3x3(3,3)
        
        Ci3x3 = (detCi3x3**(-1.0/3.0))*Ci3x3   ! C_i for the current time step
        
        Ci_new(1:2,1:2,i) = Ci3x3(1:2,1:2)
        
        call inv_matrix(C,invC)
        call inv_matrix(Cip,invCip)
        call mymulty(C_iso,invCip,multCCi)
        call dev(multCCi,devmultCCi)
        call mymulty(mu*invC,devmultCCi,Stress2PK)
        
        call mymulty(Stress2PK,trans_Fp,Couchy_tmp1)
        call mymulty(Fp,Couchy_tmp1,Couchy_tmp)   ! now Couchy_tmp is Kirchhoff
        
        do alpha=1,3
            Couchy_tmp(alpha,alpha)=Couchy_tmp(alpha,alpha)+k/10.0*(detFp**5-detFp**(-5))
        enddo
         ! now Couchy_tmp is Kirchhoff
        
        
        call inv_matrix(Fp,invFp)
        
         do alpha=1,2
            do beta=1,2
               PK1(alpha,beta,i)=0
                do gamma=1,3
                   PK1(alpha,beta,i)=PK1(alpha,beta,i)+invFp(alpha,gamma)*Couchy_tmp(gamma,beta)
                enddo
            enddo
         enddo
        
         PK1_trans_tmp(1:2,1:2)=PK1(1:2,1:2,i)   ! this way it works better :-)
         
         do alpha=1,2
          do beta=1,2
               PK1(alpha,beta,i)=PK1_trans_tmp(beta,alpha)
            enddo
        enddo
    
        Couchy(1:2,1:2,i)=Couchy_tmp(1:2,1:2)  
    enddo
return
end