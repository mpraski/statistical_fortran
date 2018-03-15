program linreg

   implicit none

   integer, parameter :: dp = selected_real_kind(15, 307)

   real(dp), allocatable :: data(:, :)

   real(dp)      :: sx = 0.0, sy = 0.0, sx2 = 0.0, sy2 = 0.0, sxy = 0.0
   real(dp)      :: m, b, r, sem, seb, sse = 0.0
   integer       :: i, n = 0, u = 10
   integer       :: status
   character(80) :: filename

   write (unit=*, fmt="(a)", advance="no") " Enter filename:  "
   read (unit=*, fmt="(a)") filename

   write (unit=*, fmt="(a)", advance="no") " Enter number of data points:  "
   read (unit=*, fmt="(i5)") n

   allocate (data(n, 2))

   open (unit=u, file=filename, action="read", status="old")

   DO i = 1, n
      read (unit=u, fmt="(2 f15.10)", iostat=status) data(i, 1), data(i, 2)
      IF (status /= 0) EXIT

      sx = sx + data(i, 1)
      sy = sy + data(i, 2)
      sx2 = sx2 + data(i, 1)*data(i, 1)
      sy2 = sy2 + data(i, 2)*data(i, 2)
      sxy = sxy + data(i, 1)*data(i, 2)
   END DO

   close (unit=u)

   m = (n*sxy - sx*sy)/(n*sx2 - sx**2)
   b = (sy*sx2 - sx*sxy)/(n*sx2 - sx**2)
   r = (sxy - sx*sy/n)/ &
       sqrt((sx2 - sx**2/n)*(sy2 - sy**2/n))

   DO i = 1, n
      sse = sse + (data(i, 2) - (b + m*data(i, 1)))**2
   END DO

   seb = sqrt(sse/(n - 2))* &
         sqrt(1/n + (sx/n)**2/ &
              (sx2 - 2*(sx/n)*sx + n*(sx/n)**2))

   sem = sqrt(sse/(n - 2))/ &
         sqrt(sx2 - 2*(sx/n)*sx + n*(sx/n)**2)

   write (unit=*, fmt="(/a,f15.10)") " Slope             m = ", m
   write (unit=*, fmt="(a, f15.10)") " y-intercept       b = ", b
   write (unit=*, fmt="(a, f15.10)") " Correlation       r = ", r
   write (unit=*, fmt="(a, f15.10)") " Standard Error of m = ", sem
   write (unit=*, fmt="(a, f15.10)") " Standard Error of b = ", seb
   write (unit=*, fmt="(a,     i5)") " Observations      n = ", n

end program linreg
