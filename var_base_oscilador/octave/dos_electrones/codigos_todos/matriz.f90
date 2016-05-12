program matriz
  implicit none
  integer, parameter :: N_dim = 49
  integer :: n1, n2, m1, m2, ind1, ind2
  real(8) :: S(N_dim**2,N_dim**2)

  S(:,:) = 0.d0;
  ind2 = 1
  do n2 = 1, N_dim
    do m2 = 1, N_dim

      ind1 = 1
      do n1 = 1, N_dim
        do m1 = 1, N_dim

          if(n2==m2) then
            if(n1==n2 .and. m1==n2) then
              S(ind1,ind2) = 1.d0;
            end if
          else
            if(n1==n2 .and. m1==m2) then
              S(ind1,ind2) = 0.5d0*sqrt(2.d0);
            elseif(n1==m2 .and. m1==n2) then
              S(ind1,ind2) = 0.5d0*sqrt(2.d0);
            end if
          end if

          ind1 = ind1 + 1;
        end do
      end do

      ind2 = ind2 + 1;
    end do
  end do

  open(11, file='matriz_cambio_base.dat')
  do ind1 = 1, N_dim**2
    write(11,6) (S(ind1,ind2), ind2 = 1, N_dim**2)
  end do
  close(11)

6 format(10000(2x,E22.14))

end program matriz
