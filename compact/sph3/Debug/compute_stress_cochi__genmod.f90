        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 11 18:02:55 2019
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
