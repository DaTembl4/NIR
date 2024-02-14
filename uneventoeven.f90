module uneven2even
implicit none
real*8, parameter, public :: a_par = 0.25, n_par = 0.5
real*8, parameter, public :: pi = acos(-1d0)
complex*16, parameter, public :: im_unit = (0, 1)
real*8, public :: h = 1d-9
integer*4, public :: m
contains 

subroutine even_split(time, magnitude, n, time_prime, magnitude_prime)
real*8, intent(in), dimension(n) :: time, magnitude
real*8, intent(out), dimension(m) :: time_prime, magnitude_prime
integer*4 :: k, j, n

magnitude_prime = 0
time_prime(1) = time(1)
magnitude_prime(1) = magnitude(1)
j = 2
do k = 2, m
    time_prime(k) = time(1) + h * (k)
    if (abs(time_prime(k) - time(j)) .lt. h) then
        magnitude_prime(k) = magnitude(j)
        j = j + 1
    endif
enddo
endsubroutine

subroutine interpolation(time, magnitude, n, time_prime, magnitude_prime, pts)
real*8, intent(in), dimension(n) :: time, magnitude
real*8, dimension(m) :: time_prime, magnitude_prime
real*8 :: step(2:n), r(3:n), s(3:n), k(2:n), l(2:n), c(2:n+1), a(2:n), b(2:n), d(2:n)
integer*4 :: n, iter, j, pts(n)

step(2:n) = time(2:n) - time(1:n-1)
s(3:n) = 2 * (step(2:n-1) + step(3:n))
r(3:n) = 3 * ((magnitude(3:n) - magnitude(2:n-1)) / step(3:n) - (magnitude(2:n-1) - magnitude(1:n-2)) / step(2:n-1))
k(2) = 0d0
l(2) = 0d0
do iter = 3, n
    l(iter) = - step(iter) / ( s(iter) + step(iter-1) * l(iter-1) )
    k(iter) = ((r(iter) - step(iter-1) * k(iter-1)) / (s(iter) + step(iter-1) * l(iter-1)))
enddo
c(2) = 0
c(n+1) = 0
do iter = 2, n - 1
    c(n-iter+1) = k(n-iter+1) - l(n-iter+1) * c(n-iter+2)
enddo
a(2:n) = magnitude(1:n-1)
d(2:n) = (c(3:n+1) - c(2:n)) / (3 * step(2:n))
b(2:n) = (magnitude(2:n) - magnitude(1:n-1)) / step(2:n) - (c(3:n+1) + 2 * c(2:n)) * step(2:n) / 3

iter = 2
do j = 1, m
    if (j .lt. pts(iter)) then
        magnitude_prime(j) = a(iter) + b(iter) * (time_prime(j) - time(iter-1)) + c(iter) * ((time_prime(j) - time(iter-1)) ** 2) &
        + d(iter) * (((time_prime(j)) - time(iter-1)) ** 3)
    else
        iter = iter + 1
    endif
enddo
endsubroutine

subroutine periodogram(n, N2, magnitude, freq, time)
real*8, intent(in), dimension(n) :: magnitude, time
real*8, allocatable :: correlogramm(:), W(:), correlogramm_prime(:)
real*8 :: sigma02, mean, freq, magnitude_centered(0:N2-1), D(0:N2/2), D_prime(0:N2/2)
integer*4 :: j, N2, n, n_star, k(0:N2-1)
complex*16 :: X(0:N2-1)

do j = 0, N2 - 1
    k(j) = j
enddo
mean = 1.0 / n * sum(magnitude)

magnitude_centered(0:n-1) = magnitude(1:n) - mean
magnitude_centered(n:N2-1) = 0

do j = 0, N2 - 1
    X(j) = sum(magnitude_centered * exp(-im_unit * 2 * pi / N2 * k * j))
enddo
do j = 0, N2 / 2
    D(j) = 1.0 / (n ** 2) * ((real(X(j))) ** 2 + (aimag(X(j))) ** 2)
enddo

sigma02 = 1.0 / (n - 1) * sum((magnitude_centered(0:n-1)) ** 2)

n_star = int(n * n_par)
allocate(correlogramm(0:n_star-1), correlogramm_prime(0:N2-1), W(0:n_star-1))

do j = 0, n_star - 1
    correlogramm(j) = 1.0 / n * real(sum(abs(X) ** 2 * exp(im_unit * 2 * pi / N2 * k * j)) / N2)
    W(j) = (1 - 2 * a_par) + 2 * a_par * cos(pi * j / n_star)
    correlogramm_prime(j) = correlogramm(j) * W(j)
enddo

correlogramm_prime(n_star:N2-1) = 0

do j = 0, N2 / 2
D_prime(j) = 1.0 / (n_star) * (2 * real(sum(correlogramm_prime * exp(-im_unit * 2 * pi / N2 * k * j)) / N2))! - correlogramm_prime(0))
enddo

print*, maxloc(D_prime) * freq / 24
print*, maxloc(D) * freq / 24

!open(3, file = 'DATA\periodogramm.dat')
open(3, file = 'periodogramm.dat')
do j = 0, N2 / 2
    write(3, *) j * freq, D(j)
enddo
!open(4, file = 'DATA\periodogramm_prime.dat')
open(4, file = 'periodogramm_prime.dat')
do j = 0, N2 / 2
    write(4, *) j * freq, D_prime(j)
enddo
!open(5, file = 'DATA\correlogramm.dat')
open(5, file = 'correlogramm.dat')
do j = 0, n_star - 1
    write(5, *) (time(j+1) - time(1)) * 24, correlogramm_prime(j)
enddo
!open(6, file = 'DATA\correlogramm_prime.dat')
open(6, file = 'correlogramm_prime.dat')
do j = 0, n_star - 1
    write(6, *) (time(j+1) - time(1)) * 24, correlogramm(j)
enddo

endsubroutine

endmodule uneven2even