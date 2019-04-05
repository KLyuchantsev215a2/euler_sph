subroutine Compute_Stress_Cochi(F,C,mu,k,N)
    
    integer:: N
    real*8 :: F(2,2,N)
    real*8 :: C(2,2,N)
    real*8 :: mu
    real*8 :: k

    real*8::detFp
    real*8 ::Fp(3,3)
    real*8 :: trans_Fp(3,3)
    real*8 ::B(3,3)
    real*8 ::B_iso(3,3)
    real*8 ::dev_B_iso(3,3)
    
    C=0
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

        do alpha=1,3
            do beta=1,3
                trans_Fp(alpha,beta)=Fp(beta,alpha)
            enddo
        enddo

        do alpha=1,3
            do beta=1,3
                B(alpha,beta)=0
                do gamma=1,3
                    B(alpha,beta)=B(alpha,beta)+Fp(alpha,gamma)*trans_Fp(gamma,beta)
                enddo
            enddo
        enddo

        B_iso=detFp**(-2.0/3.0)*B

        dev_B_iso=B_iso

        dev_B_iso(1,1)=dev_B_iso(1,1)-(1.0/3.0)*(B_iso(1,1)+B_iso(2,2)+B_iso(3,3))
        dev_B_iso(2,2)=dev_B_iso(2,2)-(1.0/3.0)*(B_iso(1,1)+B_iso(2,2)+B_iso(3,3))
        dev_B_iso(3,3)=dev_B_iso(3,3)-(1.0/3.0)*(B_iso(1,1)+B_iso(2,2)+B_iso(3,3))

        C(1:2,1:2,i)=mu*dev_B_iso(1:2,1:2)/detFp

        C(1,1,i)=C(1,1,i)+k/10.0*(detFp**5-detFp**(-5))/detFp!ÂÎÏÐÎÑ
        C(2,2,i)=C(2,2,i)+k/10.0*(detFp**5-detFp**(-5))/detFp!ÂÎÏÐÎÑ

    enddo
    
return
end