program period
implicit none
integer*4 :: i, n, m, J, jj(1)
real*8 :: delta_min, nu_max, pi = acos(-1.0), Xq = 10, g = 0.01, D_max, Dq !!! Xq рандомное, ошибка не здесь
real*8, allocatable :: t(:), f(:), h(:)
complex*16 :: im_unit = (0, 1), a
complex*16, allocatable :: D(:), W(:), C(:)

open(1, file = 'source.dat')
open(2, file = 'source_prime.dat')
read(1, *)
read(1, *) n
allocate(t(n), f(n), h(n-1))
read(1, *) t(1), f(1)
backspace(1)
do i = 1, n
    read(1, *) t(i), f(i)
    write(2, *) t(i) * 24 - floor(t(1) * 24), f(i)         !Это чтобы удобнее в гнуплоте строить график было
enddo
close(1)
close(2)

! Здесь все пункты метода CLEAN по порядку, но я не доделывал до конца, потому что возникла ошибка 

h(1:n-1) = (t(2:n) - t(1:n-1)) 

delta_min = minval(h)

nu_max = 1 / (2 * delta_min)

m = floor(4 *  nu_max * (t(n) - t(1))) !!! не знаю, почему здесь 4, это придумал не я 

allocate(D(2*m+1), W(4*m+1), C(2*m+1))

do i = 1, 2 * m + 1
    D(i) = 1.0 / n * sum(f * exp( - im_unit * 2 * pi * (((i - m - 1) / real(m)) * nu_max) * t))
enddo

do i = 1, 4 * m + 1
    W(i) = 1.0 / n * sum(exp( - im_unit * 2 * pi * (((i - 2 * m - 1) / real(m)) * nu_max) * t))
enddo

C(1:2*m+1) =  0

Dq = 1.0 / (m + 1) * sum(abs(D(m+1:2*m+1))**2) * Xq
i = 0
D_max = sqrt(maxval(abs(D)**2))

!do while (D_max .ge. Dq)
do i = 1, 3
    jj = maxloc(abs(D)**2)
    J = jj(1)
    
    D_max = sqrt(maxval(abs(D)**2))

    a = (D(J) - conjg(D(J)) * W(2*J)) / (1 - abs(W(2*J)) ** 2)

    D(1:2*m+1) = D(1:2*m+1) - g * (a * W(m+1-J:3*m+2-J) + conjg(a) * W(m+1+J:3*m+2+J))

    C(J) = C(J) + g * a
    C(-J) = C(-J) + g * conjg(a)
!    i = i + 1
    print*, D_max
enddo

end program period