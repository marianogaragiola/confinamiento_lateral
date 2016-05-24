program cb
  implicit none
  integer, parameter :: pr = kind(0.d0), N_base = 3
  integer :: ind1, ind2, n1, n2, m1, m2, N_dim
  real(pr), allocatable :: S(:,:)
  real(pr) :: kron

  N_dim = N_base**2;

  allocate(S(N_dim,N_dim))

  ind1 = 1;
  do n1 = 1, N_base
    do m1 = 1, N_base

      ind2 = 1;
      do n2 = 1, N_base
        do m2 = 1, N_base
          
          S(ind1, ind2) = sqrt(0.5_pr)*(kron(n1,n2)*kron(m1,m2) + kron(n1,m2)*kron(m1,n2))*(1._pr-kron(n2,m2)) &
                        & + kron(n1,n2)*kron(m1,n2)*kron(n2,m2);

          ind2 = ind2 + 1;
        end do
      end do

      ind1 = ind1 + 1;

    end do
  end do

  open(11, file = 'cambio_base.dat')
  do ind1 = 1, N_dim
    write(11, 6) (S(ind1, ind2), ind2 = 1, N_dim)
  end do
  close(11)

  deallocate(S)

6 format(100(2x,E14.7))

end program

function kron(n, m)
  implicit none
  integer, parameter :: pr = kind(0.d0)
  integer, intent(in) :: n, m
  real(pr) :: kron

  if(n.eq.m) then
    kron = 1._pr;
  else
    kron = 0._pr;
  end if
end function kron
