        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 26 17:32:15 2019
        MODULE PLOT_INIT__genmod
          INTERFACE 
            SUBROUTINE PLOT_INIT(X,N,COUNT_HOLE,COUNT_SECTION,          &
     &INDEX_SECTION,INDEX_HOLE)
              INTEGER(KIND=4) :: COUNT_SECTION
              INTEGER(KIND=4) :: COUNT_HOLE
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(2,N)
              INTEGER(KIND=4) :: INDEX_SECTION(COUNT_SECTION)
              INTEGER(KIND=4) :: INDEX_HOLE(COUNT_HOLE)
            END SUBROUTINE PLOT_INIT
          END INTERFACE 
        END MODULE PLOT_INIT__genmod
