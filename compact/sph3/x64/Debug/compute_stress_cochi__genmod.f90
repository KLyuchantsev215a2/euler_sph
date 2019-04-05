        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 19 09:46:58 2019
        MODULE COMPUTE_STRESS_COCHI__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_STRESS_COCHI(F,C,MU,K,N)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: F(2,2,N)
              REAL(KIND=8) :: C(2,2,N)
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: K
            END SUBROUTINE COMPUTE_STRESS_COCHI
          END INTERFACE 
        END MODULE COMPUTE_STRESS_COCHI__genmod
