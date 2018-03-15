program ols

   implicit none

   external DGELS

   integer, parameter :: dp = kind(0.d0)

   real(dp), allocatable :: data(:, :), work(:), ya(:), yab(:), xa(:, :), xab(:, :), res(:), vc(:, :), se(:), yam(:, :), vares(:, :)

   real(dp)      :: var
   integer       :: xs, xe, y, n, s, c, i, j, u = 10, status, lwork, info
   character(80) :: filename

   write (unit=*, fmt="(a)", advance="no") " Enter filename: "
   read (unit=*, fmt="(a)") filename
   write (unit=*, fmt="(a)", advance="no") " Enter number of data rows: "
   read (unit=*, fmt="(i5)") n
   write (unit=*, fmt="(a)", advance="no") " Enter number of data columns: "
   read (unit=*, fmt="(i5)") c
   write (unit=*, fmt="(a)", advance="no") " Enter index of y column: "
   read (unit=*, fmt="(i5)") y
   write (unit=*, fmt="(a)", advance="no") " Enter starting index of x column: "
   read (unit=*, fmt="(i5)") xs
   write (unit=*, fmt="(a)", advance="no") " Enter ending index of x column: "
   read (unit=*, fmt="(i5)") xe

   s = xe - xs + 2
   lwork = 2*n*s

   allocate (data(n, c), work(lwork), xa(n, s), res(n), yam(s, 1))

   open (unit=u, file=filename, action="read", status="old")

   DO i = 1, n
      read (unit=u, iostat=status, fmt=*) data(i, :)
      IF (status /= 0) EXIT
   END DO

   close (unit=u)

   ya = data(:, y)
   yab = ya
   xa(:, 1) = 1.0_dp
   xa(:, 2:s) = data(:, xs:xe)
   xab = xa

   deallocate (data)

   CALL DGELS('N', n, s, 1.0_dp, xa, n, ya, n, work, lwork, info)

   deallocate (work)

   !compute array of residual terms
   !experiment with OpenMP
   !$OMP PARALLEL DO
   DO i = 1, n
      res(i) = yab(i)
      DO j = 1, s
         res(i) = res(i) - ya(j)*xab(i, j)
      END DO
   END DO
   !$OMP END PARALLEL DO

   !variance of error terms
   var = sum(res**2)/(n - s - 1)
   !s by s variance-covariance matrix
   vc = var*inv(ata(xab))
   !standard errors of estimators
   se = sqrt(diag(vc))

   DO i = 1, s
      yam(i, 1) = ya(i)
   END DO
   !Variance of estimators
   vares = atba(yam, vc)

   write (*, *)

   DO i = 1, s
      write (unit=*, fmt="(a, i5   )", advance="no") " Estimator ", i
      write (unit=*, fmt="(a       )", advance="no") ":"
      write (unit=*, fmt="(E15.7   )", advance="no") ya(i)
      write (unit=*, fmt="(a, E15.7)", advance="no") "  (s.e. =", se(i)
      write (unit=*, fmt="(a       )") ")"
   END DO

   write (*, *)
   write (unit=*, fmt="(a, E15.7)") " Variance of error terms: ", var
   write (unit=*, fmt="(a, E15.7)") " Variance of estimators:  ", vares(1, 1)
   write (unit=*, fmt="(a, F15.7)") " R^2 (%):                 ", vares(1, 1)/(vares(1, 1) + var)*100
   write (*, *)

   IF (s <= 4) THEN
      write (unit=*, fmt="(a)") " Variance-covariance matrix of estimators: "
      DO i = 1, s
         write (*, *) (vc(i, j), j=1, s)
      END DO
   END IF

   deallocate (xa, ya, xab, yab, res, vc, se, yam)

CONTAINS
   !Get diagonals
   PURE FUNCTION diag(A) RESULT(Diags)
      real(dp), dimension(:, :), intent(in) :: A
      real(dp), dimension(size(A, 1)) :: Diags
      integer :: i

      DO i = 1, size(A, 1)
         Diags(i) = A(i, i)
      END DO
   END FUNCTION diag

   !Transpose(A) * B * A
   FUNCTION atba(A, B) RESULT(Mult)
      real(dp), dimension(:, :), intent(in) :: A, B
      real(dp), dimension(size(A, 2), size(B, 2)) :: Temp
      real(dp), dimension(size(A, 2), size(A, 2)) :: Mult
      integer :: n, m, k, j

      external DGEMM

      n = size(A, 1)
      m = size(A, 2)
      j = size(B, 1)
      k = size(B, 2)

      CALL DGEMM('T', 'N', m, k, n, 1.0_dp, A,    n, B, j, 0.0_dp, Temp, m)
      CALL DGEMM('N', 'N', m, m, k, 1.0_dp, Temp, m, A, n, 0.0_dp, Mult, m)
   END FUNCTION atba

   !Tranpose(A) * A
   FUNCTION ata(A) RESULT(Amult)
      real(dp), dimension(:, :), intent(in) :: A
      real(dp), dimension(size(A, 2), size(A, 2)) :: Amult
      integer :: n, m

      external DGEMM

      n = size(A, 1)
      m = size(A, 2)

      CALL DGEMM('T', 'N', m, m, n, 1.0_dp, A, n, A, n, 0.0_dp, Amult, m)
   END FUNCTION ata

   !Invert a matrix
   FUNCTION inv(A) RESULT(Ainv)
      real(dp), dimension(:, :), intent(in) :: A
      real(dp), dimension(size(A, 1), size(A, 2)) :: Ainv

      real(dp), dimension(size(A, 1)) :: work
      integer, dimension(size(A, 1)) :: ipiv
      integer :: n, info

      external DGETRF
      external DGETRI

      Ainv = A
      n = size(A, 1)

      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop'Matrix is singular!'
      end if

      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop'Matrix inversion failed!'
      end if
   END FUNCTION inv

end program ols
