subroutine mymulty (M,M1,multyMM1)
    
        real*8 :: M(3,3)
        real*8 :: M1(3,3)
        real*8 :: multyMM1(3,3)
        
         do alpha=1,3
            do beta=1,3
                multyMM1(alpha,beta)=0
                do gamma=1,3
                    multyMM1(alpha,beta)=multyMM1(alpha,beta)+M(alpha,gamma)*M1(gamma,beta)
                enddo
            enddo
        enddo
            
    end 