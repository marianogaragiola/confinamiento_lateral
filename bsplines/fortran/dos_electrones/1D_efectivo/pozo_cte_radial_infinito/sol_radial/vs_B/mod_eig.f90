module mod_eig
contains
subroutine eigenvalues(ndim, nev, h, ms, auval, z)
  use precision, pr => dp
  implicit none
  integer, intent(in) :: ndim, nev
  real(pr), intent(in) :: h(ndim,ndim), ms(ndim,ndim)
  real(pr), intent(out) :: auval(nev), z(ndim,nev)
  !!!!!
  integer :: ITYPE, INFO, N, LDA, LDB, m, ldz
  integer :: lwork, IL, IU
  integer :: iwork(5*ndim), ifail(ndim)
  character(1) :: JOBZ, RANGE, UPLO
  real(pr) :: vl, vu, abstol
  real(pr) :: w(ndim), work(8*ndim)


  ITYPE = 1; ! esto es A*x = (lambda)*B*x
  JOBZ = 'V'; ! solo autovalores
  RANGE = 'I';
  UPLO = 'U';
  N = ndim; LDA = ndim; LDB = ndim;
  IL = 1; IU = nev;
  abstol = 1.d-12;
  m = nev;
  ldz = ndim;
  lwork = 8*ndim;

  call dsygvx(ITYPE, JOBZ, RANGE, UPLO, N, h, LDA, ms, LDB, &
              &vl, vu, il, iu, abstol, m, w, &
              &z, ldz, work, lwork, iwork, ifail, info);

  if(0.ne.info) then
    write(*,*) " Problemas en dsygvx, info =", info
    stop
  end if

  auval(1:nev) = w(1:nev);

end subroutine eigenvalues
end module mod_eig
