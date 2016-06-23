module matrices
contains
subroutine hamiltoniano(nb, v0, me, omega, s, x_mat, p_mat, v_pozo, v_campo, ke, h, ms)
  use precision, pr => dp
  implicit none
  integer, intent(in) :: nb
  real(pr), intent(in) :: v0, me, omega
  real(pr), intent(in) :: s(nb,nb), x_mat(nb,nb), p_mat(nb,nb), v_pozo(nb,nb), v_campo(nb,nb), ke(nb,nb)
  complex(pr), intent(out) :: h(nb**3,nb**3), ms(nb**3,nb**3)
  !!!!
  integer :: i, j
  real(pr) :: Lz(nb**3,nb**3)
  real(pr) :: aux1(nb**2,nb**2), aux2(nb**3,nb**3)
  real(pr) :: T_cinetica(nb**3,nb**3), U_pozo(nb**3,nb**3), U_campo(nb**3,nb**3)

  !!!! Primero la energia cinetica
  aux1 = 0._pr; aux2 = 0._pr;
  T_cinetica = 0._pr;

  call prod_kron(nb, nb, S, S, aux1);
  call prod_kron(nb, nb*nb, ke, aux1, aux2);

  T_cinetica = T_cinetica + aux2;

  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, ke, S, aux1);
  call prod_kron(nb, nb*nb, S, aux1, aux2);

  T_cinetica = T_cinetica + aux2;

  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, S, ke, aux1);
  call prod_kron(nb, nb*nb, S, aux1, aux2);

  T_cinetica = T_cinetica + aux2;

  !!!! Segundo la energia potencial del pozo de potencial
  U_pozo = 0._pr;
  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, v_pozo, v_pozo, aux1);
  call prod_kron(nb, nb*nb, v_pozo, aux1, aux2);

  U_pozo = aux2;

  !!!! Tercero la energia potencial del campo magnetico
  U_campo = 0._pr;
  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, s, s, aux1);
  call prod_kron(nb, nb*nb, v_campo, aux1, aux2);

  U_campo = U_campo + aux2;

  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, v_campo, s, aux1);
  call prod_kron(nb, nb*nb, s, aux1, aux2);

  U_campo = U_campo + aux2;

  !!!! Cuarto la energia del termino de momento angular
  Lz = 0._pr;
  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, p_mat, s, aux1);
  call prod_kron(nb, nb*nb, x_mat, aux1, aux2);

  Lz = Lz + aux2;

  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, x_mat, s, aux1);
  call prod_kron(nb, nb*nb, p_mat, aux1, aux2);

  Lz = Lz - aux2;

  Lz = -omega*Lz;

  !!!! Por ultimo la matriz de solapamiento
  ms = 0._pr;
  aux1 = 0._pr; aux2 = 0._pr;
  call prod_kron(nb, nb, S, S, aux1);
  call prod_kron(nb, nb*nb, S, aux1, aux2);

  ms = dcmplx(aux2, 0._pr);

  h = dcmplx(T_cinetica - v0*U_pozo+ me*omega**2*U_campo, -Lz);

end subroutine hamiltoniano
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prod_kron(NA, NB, A, B, kron)
  use precision, pr => dp
  implicit none
  integer, intent(in) :: NA, NB
  real(pr), dimension(NA,NA), intent(in) ::A
  real(pr), dimensiOn(NB,NB), intent(in) :: B !dummy arguments
  real(pr), dimension(NA*NB,NA*NB) :: kron !output matrix of the kronecker product
  integer :: i,j !loop counters

  do i = 1, NA
    do j = 1, NA
      kron(1+NB*(i-1):NB+NB*(i-1),1+NB*(j-1):NB+NB*(j-1)) = A(i,j)*B
    end do
  end do

end subroutine prod_kron
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculo_matrices(kord, l_interval, n_cuad, nk, nb, me, sigma, t, k, x, w, s, x_mat, p_mat, v_pozo, v_campo, ke)
  use precision, pr => dp
  implicit none
  integer, intent(in) :: kord, l_interval, n_cuad, nk, nb
  integer, intent(in) :: k(nk)
  real(pr), intent(in) :: me, sigma
  real(pr), intent(in) :: t(nk), x(l_interval,n_cuad), w(l_interval,n_cuad)
  real(pr), intent(out) :: s(nb,nb), v_pozo(nb,nb), v_campo(nb,nb), ke(nb,nb)
  real(pr), intent(out) :: p_mat(nb,nb), x_mat(nb,nb)
  !!!
  integer :: i, j, m, n, im, in
  real(pr) :: bsp(kord)
  real(pr) :: bm, bn, zz

  s = 0._pr ; v_pozo = 0._pr; v_campo = 0._pr ; ke = 0._pr
  p_mat = 0._pr ; x_mat = 0._pr;

  do i = kord, kord+l_interval-1
    do j = 1, n_cuad
      zz = x(k(i),j);
      call bsplvb(t, kord, 1, zz, i, bsp);

      do m = 1,kord
        im = i-kord+m-1;
        if( im>0.and.im<nb+1 ) then
          do n = 1,kord
            in = i-kord+n-1;

            if( in>0.and.in<nb+1 )then
              s(im,in) = s(im,in) + bsp(m)*bsp(n)*w(k(i),j);

              v_pozo(im,in) = v_pozo(im,in) + w(k(i),j)*bsp(m)*bsp(n)*exp(-zz**2/(2._pr*sigma**2));

              v_campo(im,in) = v_campo(im,in) + 0.5_pr*w(k(i),j)*bsp(m)*bsp(n)*zz**2;

              x_mat(im,in) = x_mat(im,in) + zz*w(k(i),j)*bsp(m)*bsp(n);

            endif
          end do

          do n = i-kord+1, i
            if(n>1 .and. n<nb+2) then
              call bder(zz, t, bm, bn, kord, nk, 0, n, i)

              p_mat(im,n-1) = p_mat(im,n-1) + w(k(i),j)*bsp(m)*bn;

            end if
          end do

        endif
      end do
    end do
  end do

  do i = kord, kord+l_interval-1
    do m = i-kord+1, i
      if(m>1.and.m<nb+2)then
        do n = m, i
          if(n<nb+2)then
            do j = 1, n_cuad

              call bder(x(k(i),j), t, bm, bn, kord, nk, m, n, i)

              ke(m-1,n-1) = ke(m-1,n-1) + 0.5_pr*w(k(i),j)*bm*bn/me;

            end do
          endif
        end do
      endif
    end do
  end do

  do n = 1, nb
    do m = n, nb
      s(m,n) = s(n,m);
      v_pozo(m,n) = v_pozo(n,m);
      v_campo(m,n) = v_campo(n,m);
      x_mat(m,n) = x_mat(n,m);
      ke(m,n) = ke(n,m);
     end do
   end do

end subroutine calculo_matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  bajada de http://www.atom.physto.se/~lindroth/comp08/bget.f


       subroutine bder(rr,t,dm,dn,kord,np,indexm,indexn,left)
!     returns the value of (d/dx) Bspline(kord,index) in rr
!     The first Bspline is called spline #1 (i.e. index=1)
!     the first knot point is in t(1)
!     np= number of knot points (distinct or multiple) including
!     the ghost points: N phyical points np=N +2*(kord-1)


      implicit none
      integer i,left,kord,indexm,indexn,np
      real*8 dm,dn,t(np),Sp(kord),rr

      dm=0.d0  ; dn=0.d0
!     if rr=t(np) then the routine assumes that
!     there is kord knotpoints in the last physical point and
!     returns Bder.ne.zero if index is np-kord, or np-kord-1
      if(rr.gt.t(np)) return
      if(rr.lt.t(1)) return
!      do it=1,np
!        if(rr.ge.t(it)) left=it
!      end do

      if(abs(rr-t(np)).lt.1.d-10) then
        !if(index.lt.np-kord-1) return

        if(indexm.eq.np-kord) then
          dm=dble(kord-1)/(t(np)-t(np-kord))
        else if(indexm.eq.np-kord-1) then
          dm=-dble(kord-1)/(t(np)-t(np-kord))
        end if

       if(indexn.eq.np-kord) then
          dn=dble(kord-1)/(t(np)-t(np-kord))
        else if(indexn.eq.np-kord-1) then
          dn=-dble(kord-1)/(t(np)-t(np-kord))
        end if

        return
      end if

      call bsplvb(t,kord-1,1,rr,left,Sp)

  if(indexm-left+kord.ge.1.or.indexm-left+kord.le.kord)then
      i=indexm-left+kord
      if(i.eq.1) then
        dm=dble(kord-1)*(-Sp(i)/(t(indexm+kord)-t(indexm+1)))
      else if(i.eq.kord) then
        dm=dble(kord-1)*(Sp(i-1)/(t(indexm+kord-1)-t(indexm)))
      else
        dm=dble(kord-1)*(Sp(i-1)/(t(indexm+kord-1)-t(indexm))- &
     &  Sp(i)/(t(indexm+kord)-t(indexm+1)))
      end if
  endif
  if(indexn-left+kord.ge.1.or.indexn-left+kord.le.kord)then
      i=indexn-left+kord

      if(i.eq.1) then
        dn=dble(kord-1)*(-Sp(i)/(t(indexn+kord)-t(indexn+1)))
      else if(i.eq.kord) then
        dn=dble(kord-1)*(Sp(i-1)/(t(indexn+kord-1)-t(indexn)))
      else
        dn=dble(kord-1)*(Sp(i-1)/(t(indexn+kord-1)-t(indexn))- &
     &  Sp(i)/(t(indexn+kord)-t(indexn+1)))
      end if

  endif
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE bsplvb(t,jhigh,index,x,left,biatx)
!      INCLUDE 'standa.typ'
!      INCLUDE 'spldim.def'
!      include 'ndim_inc'
       implicit none
      integer(4),PARAMETER::JMAX=100
      integer index,jhigh,left,i,j,jp1
      real*8 t,x,biatx,deltal,deltar,saved,term
      DIMENSION biatx(jhigh),t(left+jhigh),deltal(jmax),deltar(jmax)
      SAVE deltal,deltar
      DATA j/1/

      GO TO (10,20),index
 10   j = 1
      biatx(1) = 1.d0
      IF (j .GE. jhigh) GO TO 99

 20   CONTINUE
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.d0
         DO i = 1,j
             term = biatx(i)/(deltar(i) + deltal(jp1-i))
             biatx(i) = saved + deltar(i)*term
             saved = deltal(jp1-i)*term
         END DO
         biatx(jp1) = saved
         j = jp1
         IF (j .LT. jhigh) GO TO 20
 99   RETURN
      END
end module matrices
