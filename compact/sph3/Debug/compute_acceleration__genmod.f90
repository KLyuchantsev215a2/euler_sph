        !COMPILER-GENERATED INTERFACE MODULE: Mon Apr 01 20:59:10 2019
        MODULE COMPUTE_ACCELERATION__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_ACCELERATION(N,H,DH,RHO_0,MU,K,ETA,      &
     &DAMPING,VOL,F,COUCHY,PK1,X,X_OLD,V_OLD,NABLA_W_0,NABLA_W,W,WPER1, &
     &WPER2,WPER3,WPER4,ACC,COUNT_HOLE,COUNT_SECTION,INDEX_SECTION,     &
     &INDEX_HOLE,CI,CI_NEW,TABLE)
              INTEGER(KIND=4) :: COUNT_SECTION
              INTEGER(KIND=4) :: COUNT_HOLE
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: H
              REAL(KIND=8) :: DH
              REAL(KIND=8) :: RHO_0
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: K
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: DAMPING
              REAL(KIND=8) :: VOL(N)
              REAL(KIND=8) :: F(2,2,N)
              REAL(KIND=8) :: COUCHY(2,2,N)
              REAL(KIND=8) :: PK1(2,2,N)
              REAL(KIND=8) :: X(2,N)
              REAL(KIND=8) :: X_OLD(2,N)
              REAL(KIND=8) :: V_OLD(2,N)
              REAL(KIND=8) :: NABLA_W_0(2,N,N)
              REAL(KIND=8) :: NABLA_W(2,N,N)
              REAL(KIND=8) :: W(N,N)
              REAL(KIND=8) :: WPER1(N,N)
              REAL(KIND=8) :: WPER2(N,N)
              REAL(KIND=8) :: WPER3(N,N)
              REAL(KIND=8) :: WPER4(N,N)
              REAL(KIND=8) :: ACC(2,N)
              INTEGER(KIND=4) :: INDEX_SECTION(COUNT_SECTION)
              INTEGER(KIND=4) :: INDEX_HOLE(COUNT_HOLE)
              REAL(KIND=8) :: CI(2,2,N)
              REAL(KIND=8) :: CI_NEW(3,3,N)
              REAL(KIND=8) :: TABLE(N,30)
            END SUBROUTINE COMPUTE_ACCELERATION
          END INTERFACE 
        END MODULE COMPUTE_ACCELERATION__genmod
