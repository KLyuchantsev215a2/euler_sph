subroutine trans (M,transM)
    
        real*8 :: M(3,3)
        real*8 :: transM(3,3)
        
        do alpha=1,3
            do beta=1,3
                transM(alpha,beta)=M(beta,alpha)
            enddo
        enddo
            
    end 