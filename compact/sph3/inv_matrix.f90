subroutine inv_matrix (M,invM)
        real*8:: M(3,3)
        real*8:: invM(3,3)
        
        invM(1,1)=M(2,2)/det(M)
        invM(1,2)=-M(1,2)/det(M)
        invM(2,1)=-M(2,1)/det(M)
        invM(2,2)=M(1,1)/det(M)
        invM(3,3)=1/M(3,3)  ! inverse of F ready!
    
    end 