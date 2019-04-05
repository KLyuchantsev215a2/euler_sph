        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 26 17:32:14 2019
        MODULE COMPUTE_STRESS_PK1__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_STRESS_PK1(F,C,PK1,MU,K,N)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: F(2,2,N)
              REAL(KIND=8) :: C(2,2,N)
              REAL(KIND=8) :: PK1(2,2,N)
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: K
            END SUBROUTINE COMPUTE_STRESS_PK1
          END INTERFACE 
        END MODULE COMPUTE_STRESS_PK1__genmod
