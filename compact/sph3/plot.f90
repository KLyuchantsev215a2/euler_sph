 SUBROUTINE plot(x,N,frame)
    use ogpf
    integer :: N,frame,i
    real*8 :: x(2,N,frame)
        !...............................................................................
        !Example 1: A very basic example
        !...............................................................................
        TYPE(gpf):: gp

        Real(wp):: xg(n)
        Real(wp):: yg(n)

        ! Input data
      
        
        
        call gp%animation_start(1)
        
        do i=1, frame, 600
            xg=x(1,1:N,i)
            yg=x(2,1:N,i)
            
            
            
            ! Annotation: set title, xlabel, ylabel
            CALL gp%title('Example 1')
            CALL gp%xlabel('x [mm]')
            CALL gp%ylabel('y [mm]')
            Call gp%options('set style data linespoints')
            !Call Plot to draw a vector against a vector of data
            call gp%options('set xrange[-0.1:2];set yrange [-0.1:2];')
            call gp%plot(xg, yg,'title "square coordinates" with points lt 6 lc rgb "#FF1100"')
        enddo
        
        call gp%animation_show()
    END SUBROUTINE plot