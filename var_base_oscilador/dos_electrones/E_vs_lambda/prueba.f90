program prueba
implicit none
integer :: i, N
real(8) :: x, V, dx
real(8) :: pi

pi = 2.d0*asin(1.d0);

N = 100;
dx = 20.d0/99.d0;

open(20, file = 'potencial.dat')
x = -10.d0;
do i = 1, N

  if( abs(x) < 5.d0) then
    V = sqrt(0.5d0*pi)*exp(x**2)*(1.d0-erf(abs(x)));
  else
    V = sqrt(0.5d0)/abs(x) - (sqrt(0.5d0)/abs(x))**3;
  end if

  write(20,*) x, V

  x = x + dx;

end do
close(20)

end program
