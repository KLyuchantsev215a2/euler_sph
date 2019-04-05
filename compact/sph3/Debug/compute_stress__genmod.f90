        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 05 14:30:45 2019
        MODULE COMPUTE_STRESS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_STRESS(F,C,MU,K,N)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: F(2,2,N)
              REAL(KIND=8) :: C(2,2,N)
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: K
            END SUBROUTINE COMPUTE_STRESS
          END INTERFACE 
        END MODULE COMPUTE_STRESS__genmod
