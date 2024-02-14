program main
use :: uneven2even
implicit none
real*8, allocatable :: t(:), f(:), t_p(:), f_p(:)
integer*4, allocatable :: points(:)
integer*4 :: i, n, l, Nk
real*8 :: nu

open(1, file = 'source.dat')
read(1, *)
read(1, *) n
allocate(t(n), f(n), points(n))
do i = 1, n
    read(1, *) t(i), f(i)
enddo
close(1)

m = floor((t(n) - t(1)) / h)
allocate(t_p(m), f_p(m))

call even_split(t, f, n, t_p, f_p)

l = 0
do i = 1, m
    if (f_p(i) .ne. 0) then
        l = l + 1
        points(l) = i
    endif
enddo

call interpolation(t, f, n, t_p, f_p, points)

deallocate(t, f)

Nk = m / (2 ** ceiling(log(real(n)) / log(2d0)))
n = 2 ** ceiling(log(real(n)) / log(2d0))
allocate(t(n), f(n))

do i = 1, n
    t(i) = t_p(i*Nk)
    f(i) = f_p(i*Nk)
enddo
deallocate(t_p, f_p)

!open(2, file = 'DATA\t&m_interpolated.dat')
open(2, file = 't&m_interpolated.dat')
do i = 1, n
    write(2, *) (t(i)-t(1))*24, f(i)
enddo

nu = 1d0 / (2 * (t(n) - t(1)))

call periodogram(n, n * 2, f, nu, t)

endprogram main
