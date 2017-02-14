program potencial
implicit none
integer :: N, i
real(8) :: xi, xf, x, dx, v, l, pi

pi = 2.d0*asin(1.d0);
N = 20001;
l = 10;

xi = -2000.d0;
xf = 2000.d0;

dx = (xf-xi)/real(N-1, 8);

open(11, file='potencial.dat');
do i = 1, N

  x = xi + real(i-1, 8)*dx;

  if(abs(x)<= 5.d0) then
    v = sqrt(0.5d0*pi)/l*exp(x*x)*(1.d0 - erf(abs(x)));
  else
    v = sqrt(0.5d0)/l*(1.d0/abs(x) - 0.5d0/abs(x)**3);
  end if


  write(11,*) x, v

end do

close(11)

end program
