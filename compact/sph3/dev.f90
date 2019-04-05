subroutine dev (M,devM)
            real*8:: M(3,3)
            real*8:: devM(3,3)
            
            devM=M
            
            do alpha=1,3
                devM(alpha,alpha)=devM(alpha,alpha)-(1.0/3.0)*trace(M)
            enddo
   
    end