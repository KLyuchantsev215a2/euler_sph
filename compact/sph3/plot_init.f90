 SUBROUTINE plot_init(x,N,count_hole,count_section,index_section,index_hole)
    use ogpf
    integer :: N,i,count_hole,count_section
    integer :: index_section(count_section)
    integer :: index_hole(count_hole)
    real*8 :: x(2,N)
        !...............................................................................
        !Example 1: A very basic example
        !...............................................................................
        TYPE(gpf):: gp

        Real(wp):: xg(n)
        Real(wp):: yg(n)
        Real(wp):: xh(count_hole)
        Real(wp):: yh(count_hole)
        Real(wp):: xs(count_section)
        Real(wp):: ys(count_section)

        ! Input data
      
            xg=x(1,1:N)
            yg=x(2,1:N)
            
            do i=1,count_hole
                xh(i)=x(1,index_hole(i))
                yh(i)=x(2,index_hole(i))
            enddo
            
            do i=1,count_section
                xs(i)=x(1,index_section(i))
                ys(i)=x(2,index_section(i))
            enddo
            
            ! Annotation: set title, xlabel, ylabel
            CALL gp%title('Example 1')
            CALL gp%xlabel('x [mm]')
            CALL gp%ylabel('y [mm]')
            Call gp%options('set style data linespoints')
            !Call Plot to draw a vector against a vector of data
            call gp%options('set xrange[-0.1:1.3];set yrange [-0.1:0.8];')
           ! call gp%plot(xs, ys,'title "square coordinates" with points lt 6 lc rgb "#000000"'
            call gp%plot(xh, yh,'title "square coordinates" with points lt 6 lc rgb "#000000"')
  
        
      
    END SUBROUTINE plot_init