
program logistic_classifier

   implicit none

   integer, parameter :: dp = kind(0.d0)

   real(dp) :: step = 0.1e-3_dp, threshold = 0.1e-1_dp
   real(dp), allocatable :: data(:, :), t(:), x(:, :), w(:)

   real(dp) :: fit
   integer :: n, s, it, ixs, ixe, num, i, u = 10, status, correct = 0
   character(80) :: filename

   write (unit=*, fmt="(a)", advance="no") " Enter filename: "
   read (unit=*, fmt="(a)") filename
   write (unit=*, fmt="(a)", advance="no") " Enter number of data rows: "
   read (unit=*, fmt="(i5)") n
   write (unit=*, fmt="(a)", advance="no") " Enter number of data columns: "
   read (unit=*, fmt="(i5)") s
   write (unit=*, fmt="(a)", advance="no") " Enter index of class column: "
   read (unit=*, fmt="(i5)") it
   write (unit=*, fmt="(a)", advance="no") " Enter starting index of x column: "
   read (unit=*, fmt="(i5)") ixs
   write (unit=*, fmt="(a)", advance="no") " Enter ending index of x column: "
   read (unit=*, fmt="(i5)") ixe

   num = ixe - ixs + 2

   allocate (data(n, s), x(n, num), w(num))

   open (unit=u, file=filename, action="read", status="old")

   DO i = 1, n
      read (unit=u, iostat=status, fmt=*) data(i, :)
      IF (status /= 0) EXIT
   END DO

   close (unit=u)

   t = data(:, it)
   x(:, 1) = 1.0_dp
   x(:, 2:num) = data(:, ixs:ixe)

   CALL random_number(w)

   CALL ascent(w, x, t, step, threshold)

   write (*, *)

   DO i = 1, num
      write (unit=*, fmt="(a, i5   )", advance="no") " Weight ", i
      write (unit=*, fmt="(a       )", advance="no") ":"
      write (unit=*, fmt="(E15.7   )") w(i)
   END DO

   DO i = 1, n
      fit = sum(x(i, :)*w)
      IF ((t(i) == 0.0_dp .AND. fit < 0) .OR. (t(i) == 1.0_dp .AND. fit > 0)) correct = correct + 1
   END DO

   write (*, *)

   write (unit=*, fmt="(a)", advance="no") " Accuracy: "
   write (unit=*, fmt="(f5.3)") real(correct)/real(n)

CONTAINS
   SUBROUTINE ascent(w, x, t, step, threshold)
      real(dp), dimension(:), intent(out) :: w
      real(dp), dimension(:, :), intent(in) :: x
      real(dp), dimension(:), intent(in) :: t
      real(dp), dimension(size(w)) :: d, wp
      real(dp) :: step, threshold

      wp = w

      CALL derive(w, x, t, d)

      DO WHILE (norm2(d) > threshold)
         CALL derive(wp, x, t, d)
         w = wp + step*d
         wp = w
      END DO
   END SUBROUTINE

   SUBROUTINE derive(w, x, t, d)
      real(dp), dimension(:), intent(in) :: w
      real(dp), dimension(:, :), intent(in) :: x
      real(dp), dimension(:), intent(in) :: t
      real(dp), dimension(:), intent(out) :: d
      integer :: i

      d = 0.0_dp

      !$OMP PARALLEL DO REDUCTION(+:d)
      DO i = 1, size(x, 1)
         d = d + (t(i) - sigmoid(sum(w*x(i, :))))*x(i, :)
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE

   PURE FUNCTION sigmoid(x) RESULT(r)
      real(dp), intent(in) :: x
      real(dp) :: r

      r = 1.0_dp/(1.0_dp + exp(-x))
   END FUNCTION
end program