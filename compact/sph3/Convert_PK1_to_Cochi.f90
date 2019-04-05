subroutine Convert_PK1_to_Cochi(F,Cochi,PK1,N)
    
    integer:: N
    real*8 :: F(2,2,N)
    real*8 :: Cochi(2,2,N)
    real*8 :: PK1(2,2,N)
    real*8 :: PK1_trans_tmp(2,2)
    
    real*8::detFp
    real*8 ::Fp(3,3)
    
    Cochi=0
    do i=1,N

        B=0
        Fp=0
        
        do alpha=1,2
            do beta=1,2
                Fp(alpha,beta)=F(alpha,beta,i)
            enddo
        enddo

        Fp(3,3)=1
        detFp=Fp(1,1)*Fp(2,2)-Fp(1,2)*Fp(2,1)
        
         PK1_trans_tmp(1:2,1:2)= PK1(1:2,1:2,i)
         
         do alpha=1,2
          do beta=1,2
                PK1(alpha,beta,i)=PK1_trans_tmp(beta,alpha)
            enddo
         enddo
         
        do alpha=1,2
            do beta=1,2
                Cochi(alpha,beta,i)=0
                do gamma=1,3
                    Cochi(alpha,beta,i)=Cochi(alpha,beta,i)+Fp(alpha,gamma)*PK1(gamma,beta,i)
                enddo
            enddo
        enddo
        
        Cochi(1:2,1:2,i)=PK1(1:2,1:2,i)/detFp

    enddo
    
return
end