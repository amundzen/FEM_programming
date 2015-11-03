module new_library

use precision

contains                                  
! modified by lm 11/06/01 to add use precision

  subroutine form_s(GG,ell,kappa,omega,Gamma,s)
  ! forms the s vector in bicgstab(l)
  real(iwp),intent(in)::GG(:,:),kappa; integer,intent(in)::ell
  real(iwp),intent(out)::omega,Gamma(:),s(:)
  real(iwp)::HH(ell-1,ell-1),Gamma0(ell+1),p(ell-1),q(ell-1),  &
        Gamma1(ell+1),NGamma0,NGamma1,cosine
            HH= -GG(2:ell,2:ell);call invert(HH)  
            p = matmul(HH,GG(2:ell,1))  ;   q = matmul(HH,GG(2:ell,ell+1))
            Gamma0(1) = 1.0; Gamma0(ell+1) = 0.0; Gamma0(2:ell) = p
            Gamma1(1) = 0.0; Gamma1(ell+1) = 1.0; Gamma1(2:ell) = q
            NGamma0   = dot_product(Gamma0,matmul(GG,Gamma0))
            NGamma1   = dot_product(Gamma1,matmul(GG,Gamma1))
            omega     = dot_product(Gamma0,matmul(GG,Gamma1))
            cosine = abs(omega)/sqrt(abs(NGamma0*NGamma1)); omega=omega/NGamma1
            if(cosine<kappa) omega = (kappa/cosine) * omega
            Gamma = Gamma0 - omega * Gamma1
            s(1:ell) = Gamma(2:ell+1)       ;  s(ell+1) = .0
  return
  end subroutine form_s
subroutine row_sort(a,k)
! bubble sort rows of a according to k/th column
integer,intent(inout)::a(:,:); integer::k,temporary,i=1,j,l,m,n
logical::keep_going=.true. ; m=ubound(a,1); n=ubound(a,2)
  do
    if(i==m.and..not.keep_going) exit
    keep_going=.false.
    do l=1,m-i
       if(a(l,k)>a(l+1,k)) then
        do j=1,n; temporary=a(l,j); a(l,j)=a(l+1,j);a(l+1,j)=temporary
        end do  ;  keep_going = .true.
       end if
    end do  ;  i=i+1
  end do
end subroutine row_sort
subroutine read_rest(rest,nodof,nrb,nrb1,nrb2,nrb3,nrb4,nrb5,nrb6)
implicit none
integer,intent(in):: nodof,nrb,nrb1,nrb2,nrb3,nrb4,nrb5,nrb6
integer,intent(out)::rest(:,:)
integer::i,j,l,n,if,is,it,nfd(nodof),inc 
inc= 0                         
do l=1,nrb                
  read(10,*) if,it,is,nfd   
      do i = if,it,is   ; inc = inc + 1
      rest(inc,1)=i; rest(inc,2:nodof+1)=nfd ;  end do
end do
do l=1,nrb1
  read(10,*) if,it,is,nfd   
   do n=1,11
    do i=if,it,is;inc = inc + 1  
    j=i+(n-1)*462;rest(inc,1)=j; rest(inc,2:nodof+1)=nfd; end do
   end do
end do
do l=1,nrb2
  read(10,*) if,it,is,nfd    
   do n=1,10
      do i = if,it,is;inc = inc + 1    
      j=i+(n-1)*462;rest(inc,1)=j;rest(inc,2:nodof+1)=nfd;end do
   end do
end do
do l=1,nrb3
  read(10,*) if,it,is,nfd  
   do n=1,10
      do i = if,it,is;inc = inc + 1    
      j=i+(n-1)*32;rest(inc,1)=j;rest(inc,2:nodof+1)=nfd; end do
   end do
end do
do l=1,nrb4
  read(10,*) if,it,is,nfd   
   do n=1,9
      do i = if,it,is; inc = inc + 1   
      j=i+(n-1)*32;rest(inc,1)=j;rest(inc,2:nodof+1)=nfd;end do
   end do
end do
do l=1,nrb5
  read(10,*) if,it,is,nfd   
   do n=1,9
      do i = if,it,is;inc = inc + 1   
      j=i+(n-1)*11;rest(inc,1)=j;rest(inc,2:nodof+1)=nfd;end do
   end do
end do
do l=1,nrb6
  read(10,*) if,it,is,nfd     
   do n=1,9
      do i = if,it,is;inc = inc + 1     
      j=i+(n-1)*462;rest(inc,1)=j;rest(inc,2:nodof+1)=nfd;end do
   end do
end do
return    
end subroutine read_rest 
  real function norm(x)
  !calculates the l2 norm of vector x
  real(iwp),intent(in)::x(:)
  norm=sqrt(sum(x**2))
  end function norm
 subroutine formupvw(ke,c11,c12,c21,c23,c32,c24,c42)
!   forms unsymmetrical 'stiffness' matrix in 3-d for
!   u-p-v-w version of navier stokes
 real(iwp),intent(in)::c11(:,:),c12(:,:),c21(:,:),c23(:,:),c32(:,:),     &
                  c24(:,:),c42(:,:)
 real(iwp),intent(out):: ke(:,:); integer:: nod,nodf,ntot
 nod=ubound(c11,1); nodf=ubound(c21,1); ntot=nod+nodf+nod+nod
 ke(1:nod,1:nod) = c11; ke(1:nod,nod+1:nod+nodf) = c12 
 ke(nod+1:nod+nodf,1:nod) = c21
 ke(nod+1:nod+nodf,nod+nodf+1:nod+nodf+nod) = c23
 ke(nod+1:nod+nodf,nod+nodf+nod+1:ntot) = c24 
 ke(nod+nodf+1:nod+nodf+nod,nod+1:nod+nodf) = c32
 ke(nod+nodf+1:nod+nodf+nod,nod+nodf+1:nod+nodf+nod) = c11 
 ke(nod+nodf+nod+1:ntot,nod+1:nod+nodf) = c42
 ke(nod+nodf+nod+1:ntot,nod+nodf+nod+1:ntot) = c11     
 end subroutine formupvw 


 subroutine formblock_k(bigk,km,g,g_min)
 ! forms complete n*n matrix or submatrix
 implicit none
 real(iwp),intent(in)::km(:,:); integer,intent(in):: g(:),g_min
 real(iwp),intent(out)::bigk(:,:)
 integer::i,j,ndof;       ndof=ubound(km,1)
 do i = 1 , ndof
  if(g(i)/=0) then
   do j = 1 , ndof
    if(g(j)/=0) then
     bigk(g(i)-g_min+1,g(j)-g_min+1)=bigk(g(i)-g_min+1,g(j)-g_min+1) + km(i,j)
    end if
   end do
  end if
 end do
return
end subroutine formblock_k
subroutine sparin_gauss(kv,kdiag)
! Gaussian factorisation of a skyline matrix
implicit none
real(iwp),intent(out)::kv(:) ; integer,intent(in)::kdiag(:)
real(iwp)::num,den,fac; integer::n,ii,i,j,k,l,kk,l1,l2,l3; n=ubound(kdiag,1)
 do j = 1 , n-1
    den = kv(kdiag(j))
    ii = 0 
    do i = j+1 , n
    ii = ii + 1 ;   l = kdiag(i) - ii
     if(l-kdiag(i-1)>.0) then
       num = kv(l)   ; fac = num/den ; kk = -1
       do k = i , n
          kk = kk + 1; l1=kdiag(i+kk)-kk;  l2=l1-ii; l3=kdiag(i+kk-1)
          if(l2-l3>.0) then
             kv(l1) = kv(l1) - fac*kv(l2)
          end if
       end do 
     end if
    end do
 end do
 return
end subroutine sparin_gauss
subroutine spabac_gauss(kv,loads,kdiag)
! Gaussian back-substitution on a skyline matrix
implicit none
real(iwp),intent(in)::kv(:); real(iwp),intent(inout)::loads(0:)
integer,intent(in)::kdiag(:)
real(iwp)::num,den,fac,asum;integer::i,j,l,n,ii,jj,l1,l2; n=ubound(kdiag,1)
 do j = 1 , n-1
    den = kv(kdiag(j))  ;  ii = 0
    do i = j+1 , n
       ii = ii + 1      ; l = kdiag(i) - ii
       if(l-kdiag(i-1)>.0) then
          num = kv(l)   ; fac = num/den  ; loads(i)=loads(i)-fac*loads(j)
       end if
    end do
 end do
 loads(n) = loads(n)/kv(kdiag(n))
 do i = n-1 , 1 , -1
    jj = 0     ;  asum = .0
    do j = i+1 , n
       jj = jj + 1 ; l1 = kdiag(i+jj)-jj ;  l2 = kdiag(i+jj-1)
       if(l1-l2>.0) then
          asum = asum + kv(l1) * loads(j)
       end if
    end do
    loads(i) = (loads(i) - asum)/kv(kdiag(i))
 end do
 return
end subroutine spabac_gauss  
     subroutine readbf(nf,nodof,nrb)
  !----- blocks of input of nf data  -----------------------
      implicit none
      integer,intent(out)::nf(:,:); integer,intent(in)::nodof,nrb
      integer::nfd(nodof),i,j,l,n,if,it,is,nn;nn=ubound(nf,2)
      nf=1
      do l=1,nrb
       read(10,*)if,it,is,nfd
       do i=if,min(nn,it),is 
         nf(:,i) = nf(:,i) * nfd(:) 
       end do
      end do
      n=0
      do j = 1,nn  ;    do i= 1,nodof
          if(nf(i,j) /= 0) then
           n=n+1     ;   nf(i,j)=n
          end if
      end do       ;    end do
     return
     end subroutine readbf                
subroutine vmflow(stress,dsbar,vmfl)
! Forms the von Mises flow vector ************ 2D only
implicit none
real(iwp),intent(in)::stress(:),dsbar; real(iwp),intent(out)::vmfl(:)
real(iwp)::sigm; sigm=(stress(1)+stress(2)+stress(4))/3._iwp
vmfl(1)=stress(1)-sigm;vmfl(2)=stress(2)-sigm
vmfl(3)=stress(3)*2._iwp ; vmfl(4)=stress(4)-sigm
vmfl = vmfl*1.5_iwp/dsbar
return
end subroutine vmflow
subroutine formaa(flow,rmat,modee)
! Modification to dee matrix
implicit none
real(iwp),intent(in)::flow(:),rmat(:,:); real(iwp),intent(out)::modee(:,:)
real(iwp)::flowt(1,4),flowa(4,1)   ; flowt(1,:) = flow ; flowa(:,1) = flow
modee = matmul(matmul(matmul(rmat,flowa),flowt),rmat)
return
end subroutine formaa
subroutine fmrmat(vmfl,dsbar,dlam,dee,temp,rmat)
! Forms the r matrix
implicit none
real(iwp),intent(in)::vmfl(:),dsbar,dlam,dee(:,:),temp(:,:)
real(iwp),intent(out)::rmat(:,:); real(iwp)::acat(4,4),acatc(4,4),qmat(4,4)
integer::i,j  ; real(iwp):: con
do i=1,4;do j=1,4;acat(i,j)=vmfl(i)*vmfl(j);end do;end do
acat = (temp-acat)/dsbar; acatc = matmul(dee,acat); qmat = acatc*dlam
do i=1,4; qmat(i,i)=qmat(i,i)+1.; end do
 do i=1,4
  con = qmat(i,i); qmat(i,i) = 1._iwp
  qmat(i,:) = qmat(i,:)/con
  do j=1,4
   if(j/=i) then
    con = qmat(j,i); qmat(j,i)=0.0_iwp
    qmat(j,:) = qmat(j,:) - qmat(i,:)*con
   end if
  end do
 end do
 rmat = matmul(qmat,dee)
return
end subroutine fmrmat
subroutine fmacat(vmfl,temp,acat)
! Intermediate step
implicit none
real(iwp),intent(in)::vmfl(:),temp(:,:)
real(iwp),intent(out)::acat(:,:);integer::i,j
do i=1,4; do j=1,4; acat(i,j)=vmfl(i)*vmfl(j); end do; end do
acat = temp - acat
return 
end subroutine fmacat
 subroutine formke(km,kp,c,ke,theta)
 ! ke matrix for incremental Biot
 implicit none
    real(iwp),intent(in)::km(:,:),kp(:,:),c(:,:),theta
    real(iwp),intent(out)::ke(:,:)
    integer::ndof; ndof = ubound(km,1)
    ke(:ndof, :ndof)=km ;    ke(:ndof,ndof+1:)=c
    ke(ndof+1:,:ndof)=transpose(c) ; ke(ndof+1:,ndof+1:) = -theta * kp
  return
 end subroutine formke

  subroutine bmat_nonaxi(bee,radius,coord,deriv,fun,iflag,lth)
!
!      this subroutine forms the strain-displacement matrix for
!      axisymmetric solids subjected to non-axisymmetric loading
!
 implicit none
  real(iwp),intent(in)::deriv(:,:),fun(:),coord(:,:)
  real(iwp),intent(out):: bee(:,:),radius; integer,intent(in)::iflag,lth
  integer::nod,k,l,m,n; nod = ubound(deriv , 2 ) ; bee = .0_iwp
  radius = sum(fun * coord(:,1))
    do  m=1,nod
      n=3*m      ;    k=n-1      ;    l=k-1
      bee(1,l)=deriv(1,m);   bee(2,k)=deriv(2,m);   bee(3,l)=fun(m)/radius
      bee(3,n)=iflag*lth*bee(3,l) ; bee(4,l)=deriv(2,m); bee(4,k)=deriv(1,m)
      bee(5,k)=-iflag*lth*fun(m)/radius ;   bee(5,n)=deriv(2,m)
      bee(6,l)=bee(5,k) ; bee(6,n)=deriv(1,m)-fun(m)/radius
    end do
   return
  end subroutine bmat_nonaxi

subroutine ecmat(ecm,fun,ndof,nodof)
implicit none
 real(iwp),intent(in)::fun(:); real(iwp),intent(out)::ecm(:,:)
 integer,intent(in):: nodof,ndof
 integer:: nod,i,j; real(iwp)::nt(ndof,nodof),tn(nodof,ndof)
 nod = ndof/nodof; nt = .0; tn = .0
  do i = 1 , nod ; do j = 1 , nodof
     nt((i-1)*nodof+j,j) = fun(i); tn(j,(i-1)*nodof+j) = fun(i)
  end do; end do
  ecm = matmul( nt , tn )
 return
end subroutine ecmat
subroutine vmpl(e,v,stress,pl)
! plastic matrix for a von-Mises material
implicit none
real(iwp),intent(in)::e,v,stress(:); real(iwp),intent(out)::pl(:,:)
real(iwp)::sx,sy,sz,txy,tyz,tzx,dsbar,ee,term(6)
integer::i,j,nst; nst = ubound(stress,1)
select case(nst)
 case(4)
  sx=stress(1); sy=stress(2); txy=stress(3); sz=stress(4)
  dsbar=sqrt((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+6._iwp*txy**2)/sqrt(2._iwp)
  ee=1.5_iwp*e/((1._iwp+v)*dsbar*dsbar)
  term(1)=(2._iwp*sx-sy-sz)/3._iwp; term(2)=(2._iwp*sy-sz-sx)/3._iwp
  term(3)=txy             ; term(4)=(2._iwp*sz-sx-sy)/3._iwp
 case(6)
  sx=stress(1); sy=stress(2);  sz=stress(3)
  txy=stress(4); tyz=stress(5); tzx= stress(6)
  dsbar=sqrt((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+                            &
              6._iwp*txy**2 + 6._iwp*tyz**2 + 6._iwp*tzx**2)/sqrt(2._iwp)
  ee=1.5_iwp*e/((1._iwp+v)*dsbar*dsbar)
  term(1)=(2._iwp*sx-sy-sz)/3._iwp; term(2)=(2._iwp*sy-sz-sx)/3._iwp
  term(3)=(2._iwp*sz-sx-sy)/3._iwp
  term(4) = txy;  term(5) = tyz;  term(6) = tzx
end select
  do i=1,nst;do j=1,nst;pl(i,j)=term(i)*term(j)*ee;pl(j,i)=pl(i,j);end do;end do
return
end subroutine vmpl 
 subroutine formupv(ke,c11,c12,c21,c23,c32)
!
!      this subroutine forms the unsymmetrical stiffness matrix
!      for the u-v-p version of the Navier-Stokes equations
!
 implicit none
  real(iwp),intent(in):: c11(:,:),c21(:,:),c23(:,:),c32(:,:),c12(:,:)
  real(iwp),intent(out):: ke(:,:)  ; integer::nod,nodf,ntot
   nod = ubound(c11,1); nodf = ubound(c21,1) ;ntot=nod+nodf+nod
   ke(1:nod,1:nod)=c11  ;  ke(1:nod,nod+1:nod+nodf)=c12
   ke(nod+1:nod+nodf,1:nod)=c21  ; ke(nod+1:nod+nodf,nod+nodf+1:ntot)=c23
   ke(nod+nodf+1:ntot,nod+1:nod+nodf)=c32 
   ke(nod+nodf+1:ntot,nod+nodf+1:ntot)=c11
  return
 end subroutine formupv 
 subroutine gauss_band(pb,work)
!
!      this subroutine performs gaussian reduction of an
!      unsymmetric banded matrix 'pb' .  array 'work'
!      used as working space.
!
implicit none
  real(iwp),intent(in out):: pb(:,:),work(:,:)   
  integer::n,iwp1,iq,iqp,iwp11,i,j,k,l,ip,k1  ;   real(iwp)::s
  n=ubound(pb,1); iwp1=(ubound(pb,2)-1)/2+1
  iq=2*iwp1-1  ;  iqp=iwp1 ;  iwp11=iwp1-1
   do 1 i=1,iwp11
     do 1 j=1,iq
      if(j>=iwp1+i)goto 2
       pb(i,j)=pb(i,j+iwp1-i)  
      goto 1
   2   pb(i,j)=0._iwp
      pb(n-i+1,j)=0._iwp
   1 continue
      do 3 k=1,n
      l=k+iwp1-1
      if(l>n)l=n  ;  ip=0  ;    s=1.e-10_iwp
     do i=k,l
      if(abs(pb(i,1))<=s) cycle
      s=abs(pb(i,1))    ;   ip=i
     end do
      if(ip==0)goto 5   ;    if(k==n)goto 11
      work(iwp1,k)=ip  ;     iqp=iqp-1 ;   j=iwp1+ip-k
      if(iqp<j)iqp=j    ;   if(j==iwp1)goto 6
    do  j=1,iqp
      s=pb(k,j)
      pb(k,j)=pb(ip,j)
      pb(ip,j)=s
    end do
    6 k1=k+1
      do 8 i=k1,l
      s=pb(i,1)/pb(k,1)
      do 9 j=2,iq
      if(j>iqp)goto 10
      pb(i,j-1)=pb(i,j)-s*pb(k,j)
      goto 9
   10 pb(i,j-1)=pb(i,j)
    9 continue
      pb(i,iq)=0._iwp      ;    work(i-k,k)=s
    8 continue
    3 continue
    5 write(11,'(" singular")')
   11 return
 end subroutine gauss_band
 subroutine solve_band(pb,copy,ans)
!
!      this subroutine performs the gaussian back-substitution
!      on the reduced unsymmetric band matrix 'pb'.
!
 implicit none
   real(iwp),intent(in):: pb(:,:),copy(:,:); real(iwp),intent(out)::ans(0:)
   integer::iwp1,n,n1,i,iv,l,iq,iv1     ;  real(iwp) ::s 
   iwp1= (ubound(pb,2)-1)/2+1; n=ubound(pb,1); iq=2*iwp1-1  ;   n1=n-1
   do  iv=1,n1
     i=int(copy(iwp1,iv)+.5_iwp)
     if(i/=iv) then
      s=ans(iv) ;   ans(iv)=ans(i) ;    ans(i)=s
     end if
     l=iv+iwp1-1
     if(l>n)l=n   ;    iv1=iv+1
     do  i=iv1,l
      ans(i)=ans(i)-copy(i-iv,iv)*ans(iv)
     end do
   end do
   ans(n)=ans(n)/pb(n,1)
      iv=n-1
    6 s=ans(iv)
      l=iq
      if(iv+l-1>n)l=n-iv+1
     do  i=2,l
      s=s-pb(iv,i)*ans(iv+i-1)   ;      ans(iv)=s/pb(iv,1)
     end do
      iv=iv-1
      if(iv/=0)goto 6
     return
   end subroutine solve_band                              
 subroutine formku(ku,km,g)
!
!      this subroutine assembles element matrices into symmetrical
!      global matrix(stored as an upper rectangle)
!
      real(iwp),intent(in):: km(:,:) ; real(iwp),intent(out)::ku(:,:)
      integer,intent(in):: g(:) ;integer::i,j,icd,ndof; ndof=ubound(km,1)
      do  i=1,ndof
      if(g(i)/=0) then
        do  j=1,ndof
          if(g(j)/=0) then
           icd=g(j)-g(i)+1
           if(icd>=1) ku(g(i),icd)=ku(g(i),icd)+km(i,j)
          end if
        end do
      end if
     end do
   return
 end subroutine formku
 subroutine formlump(diag,emm,g)
!
!      this subroutine forms the global mass matrix as vector diag
!
  implicit none
   real(iwp),intent(in)::emm(:,:)  ; real(iwp),intent(out)::diag(0:)
   integer,intent(in)::g(:); integer::i,ndof; ndof=ubound(emm,1)
     do  i=1,ndof;   diag(g(i))=diag(g(i))+emm(i,i) ; end do
    return
 end subroutine formlump   
   function bandwidth(g) result(nband)
   ! finds the element bandwidth from g
   implicit none     ; integer :: nband
   integer,intent(in)::g(:) 
      nband= maxval(g,1,g>0)-minval(g,1,g>0)
   end function bandwidth
   subroutine formnf(nf)
   ! reform nf
   implicit none
   integer,intent(in out)::nf(:,:)
    integer:: i,j,m
    m=0
    do j=1,ubound(nf,2)
       do i=1,ubound(nf,1)
          if(nf(i,j)/=0) then
             m=m+1; nf(i,j)=m
          end if
       end do
    end do
    return
  end subroutine formnf 
  subroutine invert(matrix)
  ! invert a small square matrix onto itself
  implicit none
  real(iwp),intent(in out)::matrix(:,:)
  integer::i,k,n; real(iwp)::con  ; n= ubound(matrix,1)
  do k=1,n
     con=matrix(k,k); matrix(k,k)=1._iwp
     matrix(k,:)=matrix(k,:)/con
     do i=1,n
        if(i/=k) then
           con=matrix(i,k); matrix(i,k)=0.0_iwp
           matrix(i,:)=matrix(i,:) - matrix(k,:)*con
        end if
     end do
  end do
  return
 end subroutine invert
 function determinant (jac) result(det)
 ! returns the determinant of a 1x1 2x2 3x3 jacobian matrix 
 implicit none    ; real(iwp) :: det
 real(iwp),intent(in)::jac(:,:); integer:: it ; it = ubound(jac,1)  
 select case (it)
   case (1)
    det=1.0_iwp
   case (2)
    det=jac(1,1)*jac(2,2) - jac(1,2) * jac(2,1)
   case (3)
    det= jac(1,1)*(jac(2,2) * jac(3,3) -jac(3,2) * jac(2,3))
    det= det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
    det= det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
   case default
    print*,' wrong dimension for jacobian matrix'
 end select
 return
 end function determinant
 subroutine deemat(dee,e,v)
 ! returns the elastic dee matrix for given ih
 ! ih=3,plane strain; =4,axisymmetry or plane strain elastoplasticity
 ! =6 , three dimensional
 implicit none
   real(iwp),intent(in)::e,v; real(iwp),intent(out)::dee(:,:)
 ! local variables
   real(iwp)::v1,v2,c,vv; integer :: i,ih;  dee=0.0  ; ih = ubound(dee,1)
         v1 = 1. - v; c = e/((1.+v)*(1.-2.*v))
   select case (ih)
          case(3)
             dee(1,1)=v1*c; dee(2,2)=v1*c; dee(1,2)=v*c; dee(2,1)=v*c
             dee(3,3)=.5_iwp*c*(1.-2.*v)
          case(4)
             dee(1,1)=v1*c; dee(2,2)=v1*c; dee(4,4)=v1*c
             dee(3,3)=.5_iwp*c*(1.-2.*v) ; dee(1,2)=v*c; dee(2,1)=v*c
             dee(1,4)=v*c; dee(4,1)=v*c; dee(2,4)=v*c; dee(4,2)=v*c
          case(6)
             v2=v/(1._iwp-v); vv=(1._iwp-2._iwp*v)/(1.-v)*.5_iwp
             do i=1,3; dee(i,i)=1._iwp;end do; do i=4,6; dee(i,i)=vv; end do
             dee(1,2)=v2; dee(2,1)=v2; dee(1,3)=v2; dee(3,1)=v2
             dee(2,3)=v2; dee(3,2)=v2
             dee = dee*e/(2._iwp*(1._iwp+v)*vv)
          case default
             print*,'wrong size for dee matrix'
   end select
 return
 end subroutine deemat    
 subroutine beemat(bee,deriv)
 ! bee matrix for 2-d elasticity or elastoplasticity (ih=3 or 4 respectively)
 ! or for 3-d (ih = 6)
 implicit none
 real(iwp),intent(in)::deriv(:,:);  real(iwp),intent(out)::bee(:,:)
 ! local variables
 integer::k,l,m,n , ih,nod; real(iwp)::x,y,z
 bee=0. ; ih = ubound(bee,1); nod = ubound(deriv,2)
     select case (ih)
       case(3,4)
        do m=1,nod
           k=2*m; l=k-1; x=deriv(1,m); y=deriv(2,m)
           bee(1,l)=x; bee(3,k)=x; bee(2,k)=y; bee(3,l)=y
        end do
       case(6)
        do m=1,nod
           n=3*m;  k=n-1; l=k-1
           x=deriv(1,m); y=deriv(2,m); z=deriv(3,m)
           bee(1,l)=x; bee(4,k)=x; bee(6,n)=x
           bee(2,k)=y; bee(4,l)=y; bee(5,n)=y
           bee(3,n)=z; bee(5,k)=z; bee(6,l)=z
        end do
       case default
        print*,'wrong dimension for nst in bee matrix'        
      end select   
  return
 end subroutine beemat
 subroutine bmataxi(bee,radius,coord,deriv,fun)
 ! b matrix for axisymmetry
 real(iwp),intent(in)::deriv(:,:),fun(:),coord(:,:)
 real(iwp),intent(out)::bee(:,:),radius
 integer::nod ,k,l,m; real(iwp) :: x,y
  radius = sum(fun * coord(:,1))  ; nod = ubound(deriv , 2) ; bee = .0_iwp
  do m = 1 , nod
   k=2*m; l = k-1 ; x = deriv(1,m); bee(1,l) = x; bee(3 , k) = x
   y = deriv(2,m); bee(2,k)=y; bee(3,l) = y; bee(4,l)=fun(m)/radius
  end do
  return
 end subroutine bmataxi

 subroutine sample(element,s,wt)
 ! returns the local coordinates of the integrating points
 implicit none 
  real(iwp),intent(out)::s(:,:),wt(:)  ; character(*),intent(in):: element
  integer::nip ;  real(iwp):: root3, r15 , w(3),v(9),b,c   
  root3 = 1./sqrt(3.)   ;  r15 = .2*sqrt(15.) 
  nip = ubound( s , 1 ) 
         w = (/5./9.,8./9.,5./9./); v=(/5./9.*w,8./9.*w,5./9.*w/)
     select case (element)
            case('line')
             select case(nip)
              case(1)
               s(1,1)=0.  ;  wt(1)=2.   
              case(2)
               s(1,1)=root3  ; s(2,1)=-s(1,1)  ;  wt(1)=1.  ; wt(2)=1.
              case(3)
               s(1,1)=r15 ; s(2,1)=.0     ; s(3,1)=-s(1,1)
               wt = w
              case(4)
               s(1,1)=.861136311594053  ; s(2,1)=.339981043584856  
               s(3,1)=-s(2,1)  ; s(4,1)=-s(1,1)
               wt(1)=.347854845137454 ; wt(2)=.652145154862546 
               wt(3)=wt(2) ; wt(4)=wt(1)
              case(5)
               s(1,1)=.906179845938664 ; s(2,1)=.538469310105683  
               s(3,1)=.0 ; s(4,1)=-s(2,1) ; s(5,1)=-s(1,1)
               wt(1)=.236926885056189 ; wt(2)=.478628670499366
               wt(3)=.568888888888889 ; wt(4)=wt(2) ; wt(5)=wt(1)
              case(6)
               s(1,1)=.932469514203152 ; s(2,1)=.661209386466265 
               s(3,1)=.238619186083197
               s(4,1)=-s(3,1) ; s(5,1)=-s(2,1) ; s(6,1)=-s(1,1)
               wt(1)=.171324492379170 ; wt(2)=.360761573048139 
               wt(3)=.467913934572691
               wt(4)=wt(3); wt(5)=wt(2) ; wt(6)=wt(1)
                    case default
                     print*,"wrong number of integrating points for a line"
             end select
            case('triangle') 
             select case(nip)
              case(1)   ! for triangles weights multiplied by .5
                s(1,1)=1./3.  ; s(1,2)=1./3.  ;  wt(1)= .5
              case(3)   
               s(1,1)=.5 ;  s(1,2)=.5 ;  s(2,1)=.5  
               s(2,2)=0.;  s(3,1)=0.  ;  s(3,2)=.5
               wt(1)=1./3.  ;  wt(2)=wt(1) ; wt(3)=wt(1)   ; wt = .5*wt
              case(6)
 s(1,1)=.816847572980459  ; s(1,2)=.091576213509771
 s(2,1)=s(1,2);  s(2,2)=s(1,1) ;  s(3,1)=s(1,2); s(3,2)=s(1,2)
 s(4,1)=.108103018168070 ;  s(4,2)=.445948490915965
 s(5,1)=s(4,2) ;   s(5,2)=s(4,1) ;  s(6,1)=s(4,2)  ; s(6,2)=s(4,2)
 wt(1)=.109951743655322 ;   wt(2)=wt(1)  ;   wt(3)=wt(1)
 wt(4)=.223381589678011 ;   wt(5)=wt(4)  ;   wt(6)=wt(4)    ; wt = .5*wt
              case(7)
 s(1,1)=1./3. ; s(1,2)=1./3.;s(2,1)=.797426985353087 ;s(2,2)=.101286507323456
 s(3,1)=s(2,2) ;  s(3,2)=s(2,1) ; s(4,1)=s(2,2) ;  s(4,2)=s(2,2)
 s(5,1)=.470142064105115 ;   s(5,2)=.059715871789770
 s(6,1)=s(5,2) ; s(6,2)=s(5,1);  s(7,1)=s(5,1);  s(7,2)=s(5,1)
 wt(1)=.225 ; wt(2)=.125939180544827 ;  wt(3)=wt(2);  wt(4)=wt(2)
 wt(5)=.132394152788506;  wt(6)=wt(5)      ;  wt(7)=wt(5)     ;wt = .5*wt
              case(12)
 s(1,1)=.873821971016996 ; s(1,2)=.063089014491502
 s(2,1)=s(1,2) ;  s(2,2)=s(1,1);  s(3,1)=s(1,2) ;  s(3,2)=s(1,2)
 s(4,1)=.501426509658179 ;  s(4,2)=.249286745170910
 s(5,1)=s(4,2); s(5,2)=s(4,1)   ;  s(6,1)=s(4,2) ;  s(6,2)=s(4,2)
 s(7,1)=.636502499121399 ;      s(7,2)=.310352451033785
 s(8,1)=s(7,1) ;  s(8,2)=.053145049844816 ;  s(9,1)=s(7,2) ; s(9,2)=s(7,1)
 s(10,1)=s(7,2) ; s(10,2)=s(8,2) ; s(11,1)=s(8,2);   s(11,2)=s(7,1)
 s(12,1)=s(8,2) ;  s(12,2)=s(7,2)
 wt(1)=.050844906370207 ; wt(2)=wt(1); wt(3)=wt(1)
 wt(4)=.116786275726379 ; wt(5)=wt(4); wt(6)=wt(4)
 wt(7)=.082851075618374 ; wt(8:12)=wt(7)           ; wt = .5*wt
              case(16)
 s(1,1)=1./3. ;  s(1,2)=1./3.  ;  s(2,1)=.658861384496478
 s(2,2)=.170569307751761 ; s(3,1)=s(2,2)   ;  s(3,2)=s(2,1)
 s(4,1)=s(2,2)  ; s(4,2)=s(2,2)
 s(5,1)=.898905543365938 ; s(5,2)=.050547228317031
 s(6,1)=s(5,2);  s(6,2)=s(5,1) ; s(7,1)=s(5,2)  ;  s(7,2)=s(5,2)
 s(8,1)=.081414823414554; s(8,2)=.459292588292723
 s(9,1)=s(8,2)  ;  s(9,2)=s(8,1);  s(10,1)=s(8,2) ;  s(10,2)=s(8,2)
 s(11,1)=.008394777409958; s(11,2)=.263112829634638
 s(12,1)=s(11,1)    ;  s(12,2)=.728492392955404
 s(13,1)=s(11,2) ;   s(13,2)=s(11,1)  ;  s(14,1)=s(11,2); s(14,2)=s(12,2)
 s(15,1)=s(12,2) ;  s(15,2)=s(11,1) ;  s(16,1)=s(12,2) ;  s(16,2)=s(11,2)
 wt(1)=.144315607677787 ; wt(2)=.103217370534718 ; wt(3)=wt(2); wt(4)=wt(2)
 wt(5)=.032458497623198 ; wt(6)=wt(5)   ;  wt(7)=wt(5)
 wt(8)=.095091634267284 ; wt(9)=wt(8)   ;  wt(10)=wt(8)
 wt(11)=.027230314174435 ; wt(12:16) = wt(11)  ;     wt = .5*wt
              case default
                  print*,"wrong number of integrating points for a triangle"
             end select
            case ('quadrilateral')
             select case (nip)
              case(1)
                s(1,1) = .0 ; wt(1) = 4.
              case(4)
                s(1,1)=-root3; s(1,2)= root3
                s(2,1)= root3; s(2,2)= root3
                s(3,1)=-root3; s(3,2)=-root3
                s(4,1)= root3; s(4,2)=-root3
                wt = 1.0
              case(9)
                s(1:7:3,1) = -r15; s(2:8:3,1) = .0
                s(3:9:3,1) =  r15; s(1:3,2)   = r15
                s(4:6,2)   =  .0 ; s(7:9,2)   =-r15
                     wt= v
              case default
                print*,"wrong number of integrating points for a quadrilateral"
            end select
            case('tetrahedron')    
             select case(nip)
              case(1)          ! for tetrahedra weights multiplied by 1/6
                 s(1,1)=.25    ; s(1,2)=.25  ;  s(1,3)=.25   ; wt(1)=1./6.
              case(4)
               s(1,1)=.58541020 ; s(1,2)=.13819660  ;  s(1,3)=s(1,2)
               s(2,2)=s(1,1) ; s(2,3)=s(1,2)  ;  s(2,1)=s(1,2)
               s(3,3)=s(1,1) ; s(3,1)=s(1,2)  ;  s(3,2)=s(1,2)
               s(4,1)=s(1,2) ; s(4,2)=s(1,2)  ;  s(4,3)=s(1,2) ; wt(1:4)=.25/6.
              case(5)
                s(1,1)=.25  ;  s(1,2)=.25   ; s(1,3)=.25 ;  s(2,1)=.5
                s(2,2)=1./6. ;  s(2,3)=s(2,2);  s(3,2)=.5
                s(3,3)=1./6.  ;   s(3,1)=s(3,3)   ;   s(4,3)=.5
                s(4,1)=1./6. ;    s(4,2)=s(4,1);    s(5,1)=1./6.
                s(5,2)=s(5,1) ;  s(5,3)=s(5,1) 
                wt(1)=-.8  ;  wt(2)=9./20. ;   wt(3:5)=wt(2)   ; wt =wt/6.  
              case(6)
         wt = 4./3.        ;  s(6,3) = 1.
         s(1,1)=-1. ;s(2,1)=1. ; s(3,2)=-1. ; s(4,2)=1. ;  s(5,3)=-1. 
              case default
               print*,"wrong number of integrating points for a tetrahedron"
            end select
            case('hexahedron')
             select case ( nip )
              case(1)
                     s(1,1) = .0 ; wt(1) = 8.    
              case(8)   
                     s(1,1)= root3;s(1,2)= root3;s(1,3)= root3
                     s(2,1)= root3;s(2,2)= root3;s(2,3)=-root3
                     s(3,1)= root3;s(3,2)=-root3;s(3,3)= root3
                     s(4,1)= root3;s(4,2)=-root3;s(4,3)=-root3
                     s(5,1)=-root3;s(5,2)= root3;s(5,3)= root3
                     s(6,1)=-root3;s(6,2)=-root3;s(6,3)= root3
                     s(7,1)=-root3;s(7,2)= root3;s(7,3)=-root3
                     s(8,1)=-root3;s(8,2)=-root3;s(8,3)=-root3
                     wt = 1.0                                                
              case(14)
          b=0.795822426     ;      c=0.758786911
          wt(1:6)=0.886426593   ; wt(7:) =  0.335180055
          s(1,1)=-b ; s(2,1)=b  ;  s(3,2)=-b ;   s(4,2)=b
          s(5,3)=-b   ;     s(6,3)=b
          s(7:,:) = c
          s(7,1)=-c  ;  s(7,2)=-c  ; s(7,3)=-c ; s(8,2)=-c ;   s(8,3)=-c
          s(9,1)=-c  ;  s(9,3)=-c  ; s(10,3)=-c; s(11,1)=-c
          s(11,2)=-c ;  s(12,2)=-c ; s(13,1)=-c
              case(15)
          b=1.     ;      c=0.674199862
          wt(1)=1.564444444 ;  wt(2:7)=0.355555556  ; wt(8:15)=0.537777778
          s(2,1)=-b  ;    s(3,1)=b  ;    s(4,2)=-b  ;    s(5,2)=b
          s(6,3)=-b  ;    s(7,3)=b  ;    s(8:,:)=c  ;    s(8,1)=-c
          s(8,2)=-c  ;    s(8,3)=-c ;    s(9,2)=-c  ;    s(9,3)=-c
          s(10,1)=-c ;    s(10,3)=-c  ;  s(11,3)=-c ;    s(12,1)=-c
          s(12,2)=-c ;    s(13,2)=-c  ;  s(14,1)=-c                          
              case(27)
                     wt = (/5./9.*v,8./9.*v,5./9.*v/)
                     s(1:7:3,1) = -r15; s(2:8:3,1) = .0
                     s(3:9:3,1) =  r15; s(1:3,3)   = r15
                     s(4:6,3)   =  .0 ; s(7:9,3)   =-r15
                     s(1:9,2)   = -r15
                     s(10:16:3,1) = -r15; s(11:17:3,1) = .0
                     s(12:18:3,1) =  r15; s(10:12,3)   = r15
                     s(13:15,3)   =  .0 ; s(16:18,3)   =-r15
                     s(10:18,2)   = .0   
                     s(19:25:3,1) = -r15; s(20:26:3,1) = .0
                     s(21:27:3,1) =  r15; s(19:21,3)   = r15
                     s(22:24,3)   =  .0 ; s(25:27,3)   =-r15
                     s(19:27,2)   =  r15
               case default
                 print*,"wrong number of integrating points for a hexahedron" 
             end select
            case default
             print*,"not a valid element type" 
     end select
   return
 end subroutine sample
 
 subroutine shape_der(der,points,i)
 implicit none
 integer,intent(in):: i; real(iwp),intent(in)::points(:,:)
 real,intent(out)::der(:,:)
  real(iwp)::eta,xi,zeta,xi0,eta0,zeta0,etam,etap,xim,xip,c1,c2,c3 ! local variables
  real(iwp):: t1,t2,t3,t4,t5,t6,t7,t8,t9 ,x2p1,x2m1,e2p1,e2m1,zetam,zetap,x,y,z
  integer :: xii(20), etai(20), zetai(20) ,l,ndim , nod   ! local variables
  ndim = ubound(der , 1); nod = ubound(der , 2)
  select case (ndim)
   case(1) ! one dimensional case
         xi=points(i,1)
     select case (nod)
         case(2)
           der(1,1)=-0.5 ; der(1,2)=0.5
         case(3)
           t1=-1.-xi ; t2=-xi  ; t3=1.-xi
           der(1,1)=-(t3+t2)/2.  ; der(1,2)=(t3+t1)    
           der(1,3)=-(t2+t1)/2.   
         case(4)
           t1=-1.-xi ; t2=-1./3.-xi ; t3=1./3.-xi ; t4=1.-xi
           der(1,1)=-(t3*t4+t2*t4+t2*t3)*9./16.     
           der(1,2)=(t3*t4+t1*t4+t1*t3)*27./16. 
           der(1,3)=-(t2*t4+t1*t4+t1*t2)*27./16. 
           der(1,4)=(t2*t3+t1*t3+t1*t2)*9./16.   
         case(5)
           t1=-1.-xi ; t2=-0.5-xi ; t3=-xi ; t4=0.5-xi ; t5=1.-xi
           der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*2./3.   
           der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*8./3.
           der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*4. 
           der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*8./3.
           der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*2./3.
     case default
       print*,"wrong number of nodes in shape_der"        
     end select
   case(2)      ! two dimensional elements
       xi=points(i,1); eta=points(i,2) ; c1=xi ; c2=eta ; c3=1.-c1-c2
       etam=.25*(1.-eta); etap=.25*(1.+eta); xim=.25*(1.-xi); xip=.25*(1.+xi)
       x2p1=2.*xi+1. ;   x2m1=2.*xi-1. ;  e2p1=2.*eta+1. ;   e2m1=2.*eta-1.
     select case (nod)
      case(3)
       der(1,1)=1.;der(1,3)=0.;der(1,2)=-1.
       der(2,1)=0.;der(2,3)=1.;der(2,2)=-1.
      case(6) 
       der(1,1)=4.*c1-1. ;  der(1,6)=4.*c2;  der(1,5)=0.  ; der(1,4)=-4.*c2
       der(1,3)=-(4.*c3-1.);  der(1,2)=4.*(c3-c1);   der(2,1)=0.
       der(2,6)=4.*c1 ; der(2,5)=4.*c2-1.; der(2,4)=4.*(c3-c2)
       der(2,3)=-(4.*c3-1.)  ; der(2,2)=-4.*c1
      case(15)                          
       t1=c1-.25  ;  t2=c1-.5 ;  t3=c1-.75   ;   t4=c2-.25
       t5=c2-.5   ;  t6=c2-.75 ;  t7=c3-.25  ;   t8=c3-.5 ;  t9=c3-.75
       der(1,1)=32./3.*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
       der(1,12)=128./3.*c2*(t2*(t1+c1)+c1*t1) ;  der(1,11)=64.*c2*t4*(t1+c1)
       der(1,10)=128./3.*c2*t4*t5  ; der(1,9)=0. ; der(1,8)=-128./3.*c2*t4*t5
       der(1,7)=-64.*c2*t4*(t7+c3) ; der(1,6)=-128./3.*c2*(t8*(t7+c3)+c3*t7)
       der(1,5)=-32./3.*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
       der(1,4)=128./3.*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
       der(1,3)=64.*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
       der(1,2)=128./3.*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
       der(1,13)=128.*c2*(c3*(t1+c1)-c1*t1) ;  der(1,15)=128.*c2*t4*(c3-c1)
       der(1,14)=128.*c2*(c3*t7-c1*(t7+c3))
       der(2,1)=0.0 ;  der(2,12)=128./3.*c1*t1*t2;  der(2,11)=64.*c1*t1*(t4+c2)
       der(2,10)=128./3.*c1*(t5*(t4+c2)+c2*t4)
       der(2,9)=32./3.*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
       der(2,8)=128./3.*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
       der(2,7)=64.*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
       der(2,6)=128./3.*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
       der(2,5)=-32./3.*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
       der(2,4)=-128./3.*c1*(t8*(t7+c3)+c3*t7)
       der(2,3)=-64.*c1*t1*(t7+c3)  ;  der(2,2)=-128./3.*c1*t1*t2
       der(2,13)=128.*c1*t1*(c3-c2)
       der(2,15)=128.*c1*(c3*(t4+c2)-c2*t4)
       der(2,14)=128.*c1*(c3*t7-c2*(c3+t7))        
      case (4)                                                              
       der(1,1)=-etam; der(1,2)=-etap; der(1,3)=etap; der(1,4)=etam
       der(2,1)=-xim; der(2,2)=xim; der(2,3)=xip; der(2,4)=-xip
      case(8)
       der(1,1)=etam*(2.*xi+eta); der(1,2)=-8.*etam*etap
       der(1,3)=etap*(2.*xi-eta); der(1,4)=-4.*etap*xi
       der(1,5)=etap*(2.*xi+eta); der(1,6)=8.*etap*etam
       der(1,7)=etam*(2.*xi-eta); der(1,8)=-4.*etam*xi
       der(2,1)=xim*(xi+2.*eta); der(2,2)=-4.*xim*eta
       der(2,3)=xim*(2.*eta-xi); der(2,4)=8.*xim*xip
       der(2,5)=xip*(xi+2.*eta); der(2,6)=-4.*xip*eta
       der(2,7)=xip*(2.*eta-xi); der(2,8)=-8.*xim*xip   
     case(9)
       etam = eta - 1.; etap = eta + 1.; xim = xi - 1.; xip = xi + 1.
       der(1,1)=.25*x2m1*eta*etam  ;   der(1,2)=-.5*x2m1*etap*etam
       der(1,3)=.25*x2m1*eta*etap  ;   der(1,4)=-xi*eta*etap
       der(1,5)=.25*x2p1*eta*etap  ;   der(1,6)=-.5*x2p1*etap*etam
       der(1,7)=.25*x2p1*eta*etam  ;   der(1,8)=-xi*eta*etam
       der(1,9)=2.*xi*etap*etam    ;   der(2,1)=.25*xi*xim*e2m1
       der(2,2)=-xi*xim*eta        ;   der(2,3)=.25*xi*xim*e2p1
       der(2,4)=-.5*xip*xim*e2p1   ;   der(2,5)=.25*xi*xip*e2p1
       der(2,6)=-xi*xip*eta        ;   der(2,7)=.25*xi*xip*e2m1
       der(2,8)=-.5*xip*xim*e2m1   ;   der(2,9)=2.*xip*xim*eta
     case default
       print*,"wrong number of nodes in shape_der"        
     end select
   case(3)  ! three dimensional elements
       xi=points(i,1); eta=points(i,2); zeta=points(i,3)
       etam=1.-eta ; xim=1.-xi;  zetam=1.-zeta
       etap=eta+1. ; xip=xi+1. ;  zetap=zeta+1.
    select case (nod)
     case(4)
      der(1:3,1:4) = .0
      der(1,1)=1.;  der(2,2)=1.  ;  der(3,3)=1.
      der(1,4)=-1. ;  der(2,4)=-1. ;  der(3,4)=-1.  
     case(8)
      der(1,1)=-.125*etam*zetam    ;   der(1,2)=-.125*etam*zetap
      der(1,3)=.125*etam*zetap     ;   der(1,4)=.125*etam*zetam
      der(1,5)=-.125*etap*zetam    ;   der(1,6)=-.125*etap*zetap
      der(1,7)=.125*etap*zetap     ;   der(1,8)=.125*etap*zetam
      der(2,1)=-.125*xim*zetam     ;   der(2,2)=-.125*xim*zetap
      der(2,3)=-.125*xip*zetap     ;   der(2,4)=-.125*xip*zetam
      der(2,5)=.125*xim*zetam      ;   der(2,6)=.125*xim*zetap
      der(2,7)=.125*xip*zetap      ;   der(2,8)=.125*xip*zetam
      der(3,1)=-.125*xim*etam      ;   der(3,2)=.125*xim*etam
      der(3,3)=.125*xip*etam       ;   der(3,4)=-.125*xip*etam
      der(3,5)=-.125*xim*etap      ;   der(3,6)=.125*xim*etap
      der(3,7)=.125*xip*etap       ;   der(3,8)=-.125*xip*etap  
 case(14) ! type 6 element
  x= points(i,1)    ;   y= points(i,2)  ;    z= points(i,3) 
  der(1,1)=((2.*x*y+2.*x*z+4.*x+y*z+y+z)*(y-1.)*(z-1.))/8.
  der(1,2)=((2.*x*y-2.*x*z-4.*x+y*z+y-z)*(y+1.)*(z-1.))/8.
  der(1,3)=((2.*x*y+2.*x*z+4.*x-y*z-y-z)*(y-1.)*(z-1.))/8.
  der(1,4)=((2.*x*y-2.*x*z-4.*x-y*z-y+z)*(y+1.)*(z-1.))/8.
  der(1,5)=-((2.*x*y-2.*x*z+4.*x-y*z+y-z)*(y-1.)*(z+1.))/8.
  der(1,6)=-((2.*x*y+2.*x*z-4.*x-y*z+y+z)*(y+1.)*(z+1.))/8.
  der(1,7)=-((2.*x*y-2.*x*z+4.*x+y*z-y+z)*(y-1.)*(z+1.))/8.
  der(1,8)=-((2.*x*y+2.*x*z-4.*x+y*z-y-z)*(y+1.)*(z+1.))/8.
  der(1,9)=-(y+1.)*(y-1.)*(z-1.)*x  ;   der(1,10)=(y+1.)*(y-1.)*(z+1.)*x
  der(1,11)=-(y-1.)*(z+1.)*(z-1.)*x ;   der(1,12)=(y+1.)*(z+1.)*(z-1.)*x
  der(1,13)=-((y+1.)*(y-1.)*(z+1.)*(z-1.))/2.
  der(1,14)=((y+1.)*(y-1.)*(z+1.)*(z-1.))/2.
  der(2,1)=((2.*x*y+x*z+x+2.*y*z+4.*y+z)*(x-1.)*(z-1.))/8.
  der(2,2)=((2.*x*y-x*z-x+2.*y*z+4.*y-z)*(x-1.)*(z-1.))/8.
  der(2,3)=((2.*x*y+x*z+x-2.*y*z-4.*y-z)*(x+1.)*(z-1.))/8.
  der(2,4)=((2.*x*y-x*z-x-2.*y*z-4.*y+z)*(x+1.)*(z-1.))/8.
  der(2,5)=-((2.*x*y-x*z+x-2.*y*z+4.*y-z)*(x-1.)*(z+1.))/8.
  der(2,6)=-((2.*x*y+x*z-x-2.*y*z+4.*y+z)*(x-1.)*(z+1.))/8.
  der(2,7)=-((2.*x*y-x*z+x+2.*y*z-4.*y+z)*(x+1.)*(z+1.))/8.
  der(2,8)=-((2.*x*y+x*z-x+2.*y*z-4.*y-z)*(x+1.)*(z+1.))/8.
  der(2,9)=-(x+1.)*(x-1.)*(z-1.)*y
  der(2,10)=(x+1.)*(x-1.)*(z+1.)*y
  der(2,11)=-((x+1.)*(x-1.)*(z+1.)*(z-1.))/2.
  der(2,12)=((x+1.)*(x-1.)*(z+1.)*(z-1.))/2.
  der(2,13)=-(x-1.)*(z+1.)*(z-1.)*y
  der(2,14)=(x+1.)*(z+1.)*(z-1.)*y
  der(3,1)=((x*y+2.*x*z+x+2.*y*z+y+4.*z)*(x-1.)*(y-1.))/8.
  der(3,2)=((x*y-2.*x*z-x+2.*y*z+y-4.*z)*(x-1.)*(y+1.))/8.
  der(3,3)=((x*y+2.*x*z+x-2.*y*z-y-4.*z)*(x+1.)*(y-1.))/8.
  der(3,4)=((x*y-2.*x*z-x-2.*y*z-y+4.*z)*(x+1.)*(y+1.))/8.
  der(3,5)=-((x*y-2.*x*z+x-2.*y*z+y-4.*z)*(x-1.)*(y-1.))/8.
  der(3,6)=-((x*y+2.*x*z-x-2.*y*z+y+4.*z)*(x-1.)*(y+1.))/8.
  der(3,7)=-((x*y-2.*x*z+x+2.*y*z-y+4.*z)*(x+1.)*(y-1.))/8.
  der(3,8)=-((x*y+2.*x*z-x+2.*y*z-y-4.*z)*(x+1.)*(y+1.))/8.
  der(3,9)=-((x+1.)*(x-1.)*(y+1.)*(y-1.))/2.
  der(3,10)=((x+1.)*(x-1.)*(y+1.)*(y-1.))/2.
  der(3,11)=-(x+1.)*(x-1.)*(y-1.)*z  ; der(3,12)=(x+1.)*(x-1.)*(y+1.)*z
  der(3,13)=-(x-1.)*(y+1.)*(y-1.)*z  ; der(3,14)=(x+1.)*(y+1.)*(y-1.)*z
    case(20)
      xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
      etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
      zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
      do l=1,20
         xi0=xi*xii(l); eta0=eta*etai(l); zeta0=zeta*zetai(l)
         if(l==4.or.l==8.or.l==16.or.l==20) then
          der(1,l)=-.5*xi*(1.+eta0)*(1.+zeta0)
          der(2,l)=.25*etai(l)*(1.-xi*xi)*(1.+zeta0)
          der(3,l)=.25*zetai(l)*(1.-xi*xi)*(1.+eta0)
         else if(l>=9.and.l<=12)then
          der(1,l)=.25*xii(l)*(1.-eta*eta)*(1.+zeta0)
          der(2,l)=-.5*eta*(1.+xi0)*(1.+zeta0)
          der(3,l)=.25*zetai(l)*(1.+xi0)*(1.-eta*eta)
         else if(l==2.or.l==6.or.l==14.or.l==18) then
          der(1,l)=.25*xii(l)*(1.+eta0)*(1.-zeta*zeta)
          der(2,l)=.25*etai(l)*(1.+xi0)*(1.-zeta*zeta)
          der(3,l)=-.5*zeta*(1.+xi0)*(1.+eta0)
         else
          der(1,l)=.125*xii(l)*(1.+eta0)*(1.+zeta0)*(2.*xi0+eta0+zeta0-1.)
          der(2,l)=.125*etai(l)*(1.+xi0)*(1.+zeta0)*(xi0+2.*eta0+zeta0-1.)
          der(3,l)=.125*zetai(l)*(1.+xi0)*(1.+eta0)*(xi0+eta0+2.*zeta0-1.)
         end if
      end do 
     case default
       print*,"wrong number of nodes in shape_der"        
   end select
  case default
   print*,"wrong number of dimensions in shape_der"
  end select
 return
 end subroutine shape_der
 
 subroutine shape_fun(fun,points,i)
  implicit none  
  integer,intent(in):: i; real(iwp),intent(in)::points(:,:)
  real(iwp),intent(out)::fun(:)
  real(iwp) :: eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3     !local variables
  real(iwp) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,x,y,z
  real(iwp) :: zeta,xi0,eta0,zeta0; integer::xii(20),etai(20),zetai(20),l,ndim,nod
        ndim = ubound(points , 2 ); nod = ubound(fun , 1 )  
    select case (ndim)
      case(1) ! one dimensional cases
           xi=points(i,1)
        select case(nod)
         case(2)
           t1=-1.-xi ; t2=1.-xi
           fun(1)=t2/2. ; fun(2)=-t1/2.
         case(3)
           t1=-1.-xi ; t2=-xi ; t3=1.-xi
           fun(1)=t2*t3/2. ; fun(2)=-t1*t3 ; fun(3)=t1*t2/2.
         case(4)
           t1=-1.-xi ; t2=-1./3.-xi ; t3=1./3.-xi ; t4=1.-xi
           fun(1)=t2*t3*t4*9./16.  ; fun(2)=-t1*t3*t4*27./16.
           fun(3)=t1*t2*t4*27./16. ; fun(4)=-t1*t2*t3*9./16.
         case(5)
           t1=-1.-xi ; t2=-0.5-xi ; t3=-xi ; t4=0.5-xi ; t5=1.-xi
           fun(1)=t2*t3*t4*t5*2./3. ; fun(2)=-t1*t3*t4*t5*8./3.
           fun(3)=t1*t2*t4*t5*4. ; fun(4)=-t1*t2*t3*t5*8./3.
           fun(5)=t1*t2*t3*t4*2./3.
          case default
             print*,"wrong number of nodes in shape_fun"
        end select
      case(2) ! two dimensional cases
           c1=points(i,1); c2=points(i,2); c3=1.-c1-c2 
           xi=points(i,1);  eta=points(i,2)
           etam=.25*(1.-eta); etap=.25*(1.+eta)
           xim=.25*(1.-xi); xip=.25*(1.+xi)
        select case(nod)
          case(3)
            fun = (/c1,c3,c2/)  
          case(6)
            fun(1)=(2.*c1-1.)*c1 ;  fun(6)=4.*c1*c2 ;  fun(5)=(2.*c2-1.)*c2
            fun(4)=4.*c2*c3      ;  fun(3)=(2.*c3-1.)*c3 ; fun(2)=4.*c3*c1
          case(15)
            t1=c1-.25  ;  t2=c1-.5 ;  t3=c1-.75   ;   t4=c2-.25
            t5=c2-.5   ;  t6=c2-.75 ;  t7=c3-.25  ;   t8=c3-.5 ;  t9=c3-.75
            fun(1)=32./3.*c1*t1*t2*t3   ;  fun(12)=128./3.*c1*c2*t1*t2
            fun(11)=64.*c1*c2*t1*t4     ;  fun(10)=128./3.*c1*c2*t4*t5
            fun(9)=32./3.*c2*t4*t5*t6   ;  fun(8)=128./3.*c2*c3*t4*t5
            fun(7)=64.*c2*c3*t4*t7      ;  fun(6)=128./3.*c2*c3*t7*t8
            fun(5)=32./3.*c3*t7*t8*t9   ;  fun(4)=128./3.*c3*c1*t7*t8
            fun(3)=64.*c3*c1*t1*t7      ;  fun(2)=128./3.*c3*c1*t1*t2
            fun(13)=128.*c1*c2*t1*c3    ;  fun(15)=128.*c1*c2*c3*t4
            fun(14)=128.*c1*c2*c3*t7      
          case(4)
            fun=(/4.*xim*etam,4.*xim*etap,4.*xip*etap,4.*xip*etam/)
          case(8)
            fun=(/4.*etam*xim*(-xi-eta-1.),32.*etam*xim*etap,&
                  4.*etap*xim*(-xi+eta-1.),32.*xim*xip*etap, &
                  4.*etap*xip*(xi+eta-1.), 32.*etap*xip*etam,&
                  4.*xip*etam*(xi-eta-1.), 32.*xim*xip*etam/)
          case(9)
            etam = eta - 1.; etap= eta + 1.; xim = xi - 1.; xip = xi + 1.
            fun=(/.25*xi*xim*eta*etam,-.5*xi*xim*etap*etam,&
                  .25*xi*xim*eta*etap,-.5*xip*xim*eta*etap,&
                  .25*xi*xip*eta*etap,-.5*xi*xip*etap*etam,&
                  .25*xi*xip*eta*etam,-.5*xip*xim*eta*etam,xip*xim*etap*etam/)
          case default
             print*,"wrong number of nodes in shape_fun"
        end select
      case(3) ! three dimensional cases
       xi=points(i,1); eta=points(i,2); zeta=points(i,3)
       etam=1.-eta ;  xim=1.-xi  ;  zetam=1.-zeta
       etap=eta+1. ;  xip=xi+1.   ;  zetap=zeta+1.
       select case(nod)
        case(4)
         fun(1)=xi   ;   fun(2)= eta ;  fun(3)=zeta 
         fun(4)=1.-fun(1)-fun(2)-fun(3)
        case(8)
         fun=(/.125*xim*etam*zetam,.125*xim*etam*zetap,.125*xip*etam*zetap,&
               .125*xip*etam*zetam,.125*xim*etap*zetam,.125*xim*etap*zetap,&
               .125*xip*etap*zetap,.125*xip*etap*zetam/)
        case(14) !type 6 element
    x = points(i,1);  y = points(i,2);  z = points(i,3)
  fun(1)=((x*y+x*z+2.*x+y*z+2.*y+2.*z+2.)*(x-1.)*(y-1.)*(z-1.))/8.
  fun(2)=((x*y-x*z-2.*x+y*z+2.*y-2.*z-2.)*(x-1.)*(y+1.)*(z-1.))/8.
  fun(3)=((x*y+x*z+2.*x-y*z-2.*y-2.*z-2.)*(x+1.)*(y-1.)*(z-1.))/8.
  fun(4)=((x*y-x*z-2.*x-y*z-2.*y+2.*z+2.)*(x+1.)*(y+1.)*(z-1.))/8.
  fun(5)=-((x*y-x*z+2.*x-y*z+2.*y-2.*z+2.)*(x-1.)*(y-1.)*(z+1.))/8.
  fun(6)=-((x*y+x*z-2.*x-y*z+2.*y+2.*z-2.)*(x-1.)*(y+1.)*(z+1.))/8.
  fun(7)=-((x*y-x*z+2.*x+y*z-2.*y+2.*z-2.)*(x+1.)*(y-1.)*(z+1.))/8.
  fun(8)=-((x*y+x*z-2.*x+y*z-2.*y-2.*z+2.)*(x+1.)*(y+1.)*(z+1.))/8.
  fun(9)=-((x+1.)*(x-1.)*(y+1.)*(y-1.)*(z-1.))/2.
  fun(10)=((x+1.)*(x-1.)*(y+1.)*(y-1.)*(z+1.))/2.
  fun(11)=-((x+1.)*(x-1.)*(y-1.)*(z+1.)*(z-1.))/2.
  fun(12)=((x+1.)*(x-1.)*(y+1.)*(z+1.)*(z-1.))/2.
  fun(13)=-((x-1.)*(y+1.)*(y-1.)*(z+1.)*(z-1.))/2.
  fun(14)=((x+1.)*(y+1.)*(y-1.)*(z+1.)*(z-1.))/2.      
        case(20)
           xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
           etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
           zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
           do l=1,20
            xi0=xi*xii(l); eta0=eta*etai(l); zeta0=zeta*zetai(l)
            if(l==4.or.l==8.or.l==16.or.l==20) then
              fun(l)=.25*(1.-xi*xi)*(1.+eta0)*(1.+zeta0)
            else if(l>=9.and.l<=12)then
              fun(l)=.25*(1.+xi0)*(1.-eta*eta)*(1.+zeta0)
            else if(l==2.or.l==6.or.l==14.or.l==18) then
              fun(l)=.25*(1.+xi0)*(1.+eta0)*(1.-zeta*zeta)
            else
              fun(l)=.125*(1.+xi0)*(1.+eta0)*(1.+zeta0)*(xi0+eta0+zeta0-2)
            end if
           end do
          case default
           print*,"wrong number of nodes in shape_fun"
        end select
      case default
        print*,"wrong number of dimensions in shape_fun"  
    end select
  return
 end subroutine shape_fun 
subroutine rod_km(km,ea,length)
  implicit none
  real(iwp),intent(in):: ea,length; real(iwp),intent(out)::km(:,:)
  km(1,1)=1._iwp ; km(2,2)=1._iwp; km(1,2)=-1._iwp ; km(2,1)=-1._iwp
  km=km*ea/length
end subroutine rod_km 
subroutine axialkm(km,e,area,length)
  implicit none
  real(iwp),intent(in):: e,area,length; real(iwp),intent(out)::km(:,:)
  km(1,1)=1._iwp ; km(2,2)=1._iwp; km(1,2)=-1._iwp ; km(2,1)=-1._iwp
  km=km*e*area/length
end subroutine axialkm    
  subroutine mocouf(phi,c,sigm,dsbar,theta,f)
!
!      this subroutine calculates the value of the yield function
!      for a Mohr-Coulomb material (phi in degrees)
!
 implicit none
  real(iwp),intent(in)::phi,c,sigm,dsbar,theta   ; real(iwp),intent(out)::f
  real(iwp)::phir,snph,csph,csth,snth
  phir=phi*4._iwp*atan(1.)/180._iwp
  snph=sin(phir) ;  csph=cos(phir) ; csth=cos(theta);  snth=sin(theta)
  f=snph*sigm+dsbar*(csth/sqrt(3._iwp)-snth*snph/3._iwp)-c*csph
  return
  end  subroutine mocouf
  subroutine mocouq(psi,dsbar,theta,dq1,dq2,dq3)
!
!      this subroutine forms the derivatives of a Mohr-Coulomb
!      potential function with respect to the three invariants
!      psi in degrees
!
 implicit none
  real(iwp),intent(in)::psi,dsbar,theta; real(iwp),intent(out)::dq1,dq2,dq3
  real(iwp)::psir,snth,snps,sq3,c1,csth,cs3th,tn3th,tnth     
  psir=psi*4._iwp*atan(1._iwp)/180._iwp ;  snth=sin(theta) ;   snps=sin(psir)
  sq3=sqrt(3._iwp)  ; dq1=snps
      if(abs(snth).gt..49_iwp)then
        c1=1._iwp
        if(snth.lt.0._iwp)c1=-1._iwp
        dq2=(sq3*.5_iwp-c1*snps*.5_iwp/sq3)*sq3*.5_iwp/dsbar ;   dq3=0._iwp
      else
        csth=cos(theta); cs3th=cos(3.*theta);tn3th=tan(3.*theta); tnth=snth/csth
        dq2=sq3*csth/dsbar*((1._iwp+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*.5_iwp
        dq3=1.5_iwp*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar)
      end if
  return
  end subroutine mocouq
  subroutine mocopl(phi,psi,e,v,stress,pl)
!
!      this subroutine forms the plastic stress/strain matrix
!      for a mohr-coulomb material  (phi,psi in degrees)
!
  implicit none
   real(iwp),intent(in)::stress(:),phi,psi,e,v ;real(iwp),intent(out)::pl(:,:)
   integer::i,j;  real(iwp)::row(4),col(4),sx,sy,txy,sz,pi,phir,psir,&
                            dx,dy,dz,d2,d3,th,snth,sig,rph,rps,cps,snps,sq3,&
                            cc,cph ,alp,ca,sa,dd ,snph,ee,s1,s2
   sx=stress(1);  sy=stress(2);  txy=stress(3) ;  sz=stress(4)
   pi=4._iwp*atan(1._iwp) ; phir=phi*pi/180._iwp; psir=psi*pi/180._iwp
   snph=sin(phir)  ; snps=sin(psir) ; sq3=sqrt(3._iwp) ; cc=1._iwp-2._iwp*v
   dx=(2.*sx-sy-sz)/3.    ;   dy=(2.*sy-sz-sx)/3. ;  dz=(2.*sz-sx-sy)/3.
   d2=sqrt(-dx*dy-dy*dz-dz*dx+txy*txy) ;  d3=dx*dy*dz-dz*txy*txy
   th=-3._iwp*sq3*d3/(2._iwp*d2**3)
   if(th.gt.1.)th=1.   ;  if(th.lt.-1.)th=-1.
   th=asin(th)/3.   ;    snth=sin(th)
   if(abs(snth).gt..49)then
      sig=-1.
      if(snth.lt.0.)sig=1.
      rph=snph*(1.+v)/3.;  rps=snps*(1.+v)/3. ; cps=.25*sq3/d2*(1.+sig*snps/3.)
      cph=.25*sq3/d2*(1.+sig*snph/3.)
      col(1)=rph+cph*((1.-v)*dx+v*(dy+dz));col(2)=rph+cph*((1.-v)*dy+v*(dz+dx))
      col(3)=cph*cc*txy   ; col(4)=rph+cph*((1.-v)*dz+v*(dx+dy))
      row(1)=rps+cps*((1.-v)*dx+v*(dy+dz));row(2)=rps+cps*((1.-v)*dy+v*(dz+dx))
      row(3)=cps*cc*txy ;  row(4)=rps+cps*((1.-v)*dz+v*(dx+dy))
      ee=e/((1.+v)*cc*(rph*snps+2.*cph*cps*d2*d2*cc))
  else
      alp=atan(abs((sx-sy)/(2.*txy))) ; ca=cos(alp);  sa=sin(alp)
      dd=cc*sa ;  s1=1.  ;  s2=1.
      if((sx-sy).lt..0)s1=-1.
      if(txy.lt..0)s2=-1.
      col(1)=snph+s1*dd ; col(2)=snph-s1*dd;col(3)=s2*cc*ca; col(4)=2.*v*snph
      row(1)=snps+s1*dd ;row(2)=snps-s1*dd; row(3)=s2*cc*ca ; row(4)=2.*v*snps
      ee=e/(2.*(1.+v)*cc*(snph*snps+cc))
  end if
  do  i=1,4; do j=1,4 ; pl(i,j)=ee*row(i)*col(j) ; end do ; end do
  return
 end subroutine mocopl
subroutine formkv(bk,km,g,n)
!global stiffness matrix stored as a vector (upper triangle)
implicit none
real(iwp),intent(in)::km(:,:);real(iwp),intent(out)::bk(:)
integer,intent(in)::g(:),n
integer::idof,i,j,icd,ival
idof=size(km,1)
     do i=1,idof
        if(g(i)/=0) then
           do j=1,idof
              if(g(j)/=0) then
                 icd=g(j)-g(i)+1
                 if(icd-1>=0) then
                    ival=n*(icd-1)+g(i)      
                    bk(ival)=bk(ival)+km(i,j)
                 end if
               end if
            end do
         end if
     end do
return
end subroutine formkv
subroutine fsparv(bk,km,g,kdiag)
! assembly of element matrices into skyline global matrix
implicit none
real(iwp),intent(in)::km(:,:); integer,intent(in)::g(:),kdiag(:)
real(iwp),intent(out)::bk(:) ;  integer::i,idof,k,j,iw,ival
 idof=ubound(g,1)
   do i=1,idof
      k=g(i)
      if(k/=0) then
         do j=1,idof
            if(g(j)/=0) then
               iw=k-g(j)
               if(iw>=0) then
                   ival=kdiag(k)-iw
                   bk(ival)=bk(ival)+km(i,j) 
                end if
            end if
         end do
      end if
   end do
 return
end subroutine fsparv
subroutine banred(bk,n)
! gaussian reduction on a vector stored as an upper triangle
implicit none
real(iwp),intent(in out)::bk(:);integer,intent(in)::n
integer::i,il1,kbl,j,ij,nkb,m,ni,nj,iw ; real(iwp)::sum
 iw = ubound(bk,1)/n-1
       do i=2,n
          il1=i-1;kbl=il1+iw+1
          if(kbl-n>0)kbl=n
          do j=i,kbl
             ij=(j-i)*n+i;sum=bk(ij);nkb=j-iw
             if(nkb<=0)nkb=1
             if(nkb-il1<=0)then
                do m=nkb,il1
                   ni=(i-m)*n+m ; nj=(j-m)*n+m
                   sum=sum-bk(ni)*bk(nj)/bk(m) 
                end do
             end if
             bk(ij)=sum
           end do
       end do
return
end subroutine banred   
subroutine bacsub(bk,loads)
! performs the complete gaussian backsubstitution
implicit none
real(iwp),intent(in)::bk(:);real(iwp),intent(in out)::loads(0:)
integer::nkb,k,i,jn,jj,i1,n,iw;real(iwp)::sum
n = ubound(loads,1); iw = ubound(bk,1)/n - 1
loads(1)=loads(1)/bk(1)
   do i=2,n
      sum=loads(i);i1=i-1 ; nkb=i-iw
      if(nkb<=0)nkb=1
      do k=nkb,i1
         jn=(i-k)*n+k;sum=sum-bk(jn)*loads(k)
      end do
      loads(i)=sum/bk(i)
   end do
   do jj=2,n
      i=n-jj+1;sum=.0;i1=i+1;nkb=i+iw
      if(nkb-n>0)nkb=n
      do k=i1,nkb
           jn=(k-i)*n+i  ; sum=sum+bk(jn)*loads(k)
      end do
      loads(i)=loads(i)-sum/bk(i)
   end do
return
end subroutine bacsub                                                       
subroutine print_vector(vec,n,channel)
! prints out first n terms of a vector in e format to channel x 
implicit none
real(iwp),intent(in out)::vec(0:);integer,intent(in)::n,channel
integer::i
write(channel,'(1x,6g12.4)')(vec(i),i=1,n)
return
end subroutine print_vector 
subroutine print_array(array,m,n,channel)
! prints out first m rows and n columns of 'array' to channel x
implicit none
real(iwp),intent(in out)::array(:,:);integer,intent(in)::m,n,channel
integer::i,j
   do i=1,m
      write(channel,'(1x,6g12.4)')(array(i,j),j=1,n)
   end do
end subroutine print_array  
subroutine cholin(kb)
! Choleski reduction on kb(l,iw+1) stored as a lower triangle
implicit none
real(iwp),intent(in out)::kb(:,:);integer::i,j,k,l,ia,ib,n,iw;real(iwp)::x
n=ubound(kb,1); iw=ubound(kb,2)-1
   do i=1,n
      x=.0
      do j=1,iw; x=x+kb(i,j)**2; end do
      kb(i,iw+1)=sqrt(kb(i,iw+1)-x)
      do k=1,iw
         x=.0
         if(i+k<=n) then
           if(k/=iw) then
             do l=iw-k,1,-1
                x=x+kb(i+k,l)*kb(i,l+k)
             end do
           end if
           ia=i+k; ib=iw-k+1
           kb(ia,ib)=(kb(ia,ib)-x)/kb(i,iw+1)
         end if
      end do
   end do
 return
end subroutine cholin  
subroutine chobac(kb,loads)
!Choleski back-substitution
implicit none
real(iwp),intent(in)::kb(:,:);real(iwp),intent(in out)::loads(0:)
integer::iw,n,i,j,k,l,m; real(iwp)::x
n=size(kb,1); iw=size(kb,2)-1
loads(1)=loads(1)/kb(1,iw+1)
    do i=2,n
       x=.0;k=1
       if(i<=iw+1)k=iw-i+2
       do j=k,iw; x=x+kb(i,j)*loads(i+j-iw-1); end do
       loads(i)=(loads(i)-x)/kb(i,iw+1)
    end do
    loads(n)=loads(n)/kb(n,iw+1)
    do i=n-1,1,-1
       x=0.0; l=i+iw
       if(i>n-iw)l=n;  m=i+1
       do j=m,l; x=x+kb(j,iw+i-j+1)*loads(j); end do
       loads(i)=(loads(i)-x)/kb(i,iw+1)
    end do
  return
end subroutine chobac
subroutine sparin(a,kdiag)
! Choleski factorisation of variable bandwidth matrix a
! stored as a vector and overwritten
implicit none
real(iwp),intent(in out)::a(:);integer,intent(in)::kdiag(:)
integer::n,i,ki,l,kj,j,ll,m,k; real(iwp)::x
 n=ubound(kdiag,1)  ; a(1)=sqrt(a(1))
 do i=2,n
    ki=kdiag(i)-i;  l=kdiag(i-1)-ki+1
    do j=l,i
       x=a(ki+j);  kj=kdiag(j)-j
       if(j/=1) then
          ll=kdiag(j-1)-kj+1; ll=max0(l,ll)
          if(ll/=j) then
              m=j-1
              do k=ll,m ; x=x-a(ki+k)*a(kj+k) ; end do
          end if
       end if
       a(ki+j)=x/a(kj+j)
    end do
    a(ki+i)=sqrt(x)
 end do
 return
end subroutine sparin
subroutine spabac(a,b,kdiag)
! Choleski forward and backward substitution combined
! variable bandwidth factorised matrix a stored as a vector 
implicit none
real(iwp),intent(in)::a(:);real(iwp),intent(in out)::b(0:)
integer,intent(in)::kdiag(:)
integer::n,i,ki,l,m,j,it,k; real(iwp)::x
n=ubound(kdiag,1)
 b(1)=b(1)/a(1)
  do i=2,n
     ki=kdiag(i)-i;  l=kdiag(i-1)-ki+1 ; x=b(i)
     if(l/=i) then
        m=i-1
        do j=l,m ; x=x-a(ki+j)*b(j); end do
     end if
     b(i)=x/a(ki+i)
  end do
  do it=2,n
     i=n+2-it; ki=kdiag(i)-i; x=b(i)/a(ki+i); b(i)=x; l=kdiag(i-1)-ki+1
     if(l/=i) then
       m=i-1
       do k=l,m; b(k)=b(k)-x*a(ki+k); end do
     end if
  end do
 b(1)=b(1)/a(1)
 return
end subroutine spabac               
subroutine checon(loads,oldlds,tol,converged)
! sets converged to .false. if relative change in loads and
! oldlds is greater than tol and updates oldlds
implicit none
real(iwp),intent(in)::loads(0:),tol;real(iwp),intent(in out)::oldlds(0:)
logical,intent(out)::converged
  converged=.true.
  converged=(maxval(abs(loads-oldlds))/maxval(abs(loads))<=tol)
  oldlds=loads
 return
end subroutine checon
subroutine formkb(kb,km,g)
! lower triangular global stiffness kb stored as kb(n,iw+1)
implicit none
real(iwp),intent(in)::km(:,:);real(iwp),intent(out)::kb(:,:)
integer,intent(in)::g(:);integer::iw,idof,i,j,icd
idof=size(km,1);  iw=size(kb,2)-1
   do i=1,idof
      if(g(i)>0) then 
         do j=1,idof
            if(g(j)>0) then
               icd=g(j)-g(i)+iw+1
               if(icd-iw-1<=0) kb(g(i),icd)= kb(g(i),icd) +km(i,j)
            end if
         end do
      end if
   end do
 return
end subroutine formkb
subroutine formtb(kb,km,g)
! assembles unsymmetrical band matrix kb from constituent km
  implicit none
  real(iwp),intent(in)::km(:,:); integer,intent(in)::g(:)
  real(iwp),intent(out)::kb(:,:); integer::i,j,idof,icd,iw
  idof=size(km,1); iw=(size(kb,2)-1)/2
  do i=1,idof
     if(g(i)/=0) then
        do j=1,idof
           if(g(j)/=0) then
              icd=g(j)-g(i)+iw+1
              kb(g(i),icd)=kb(g(i),icd)+km(i,j)
           end if
        end do
     end if
  end do
 return
end subroutine formtb
subroutine bantmul(kb,loads,ans)
! multiplies unsymmetrical band kb by vector loads
! could be much improved for vector processors
! look out for chance for zero-sized arrays
implicit none
  real(iwp),intent(in)::kb(:,:),loads(0:); real(iwp),intent(out)::ans(0:)
  integer::i,j,k,l,m,n,iw; real(iwp)::x
  n=size(kb,1); l=size(kb,2); iw=(l-1)/2
    do i=1,n
       x=.0; k=iw+2
       do j=1,l
          k=k-1; m=i-k+1
          if(m<=n.and.m>=1)x=x+kb(i,j)*loads(m)
       end do
       ans(i)=x
    end do
 return
end subroutine bantmul 
subroutine linmul(bk,disps,loads)
! matrix-vector multiply for symmetric matrix bk
! stored in upper triangular form as a vector
implicit none
  real(iwp),intent(in)::bk(:),disps(0:); real(iwp),intent(out)::loads(0:)
   integer::i,j,n,iw; real(iwp)::x; n=ubound(disps,1);iw=ubound(bk,1)/n-1 
    do i=1,n
       x=0.0
       do j=1,iw+1
          if(i+j<=n+1)    x=x+bk(n*(j-1)+i)*disps(i+j-1)
       end do
       do j=2,iw+1
          if(i-j+1>=1)    x=x+bk((n-1)*(j-1)+i)*disps(i-j+1)
       end do    
       loads(i)=x
    end do
  return
end subroutine linmul
subroutine linmul_sky(bp,disps,loads,kdiag)
! skyline product of symmetric matrix and a vector
implicit none
 real(iwp),intent(in)::bp(:),disps(0:);real(iwp),intent(out)::loads(0:)
 integer,intent(in)::kdiag(:); integer::n,i,j,low,lup,k; real(iwp)::x
  n=ubound(disps,1)
  do i = 1 , n
     x = .0 ; lup=kdiag(i)
     if(i==1)low=lup; if(i/=1)low=kdiag(i-1)+1
     do j = low , lup
        x = x + bp(j) * disps(i + j - lup) 
     end do
     loads(i) = x
     if(i == 1) cycle   ; lup = lup - 1
     do j = low , lup
        k = i + j -lup - 1
        loads(k) = loads(k) + bp(j)*disps(i)        
     end do
  end do
 return
end subroutine linmul_sky  
function formxi(fsoil,fmax,rf,rm,ro) result(xi)
! soil spring stiffness in coupled pile-soil analysis
implicit none     ; real(iwp) :: xi
 real(iwp),intent(in)::fsoil,fmax,rf,rm,ro;  real(iwp)::phi
 phi = fsoil*ro*rf/fmax
 xi = log((rm-phi)/(ro-phi))+phi*(rm-ro)/((rm-phi)*(ro-phi))
return
end function formxi  
subroutine fkdiag(kdiag,g)
! finds the maximum bandwidth for each freedom
implicit none
 integer,intent(in)::g(:); integer,intent(out)::kdiag(:)
 integer::idof,i,iwp1,j,im,k
  idof=size(g)
  do i = 1,idof
     iwp1=1
     if(g(i)/=0) then
        do j=1,idof
           if(g(j)/=0) then
              im=g(i)-g(j)+1
              if(im>iwp1) iwp1=im
           end if
        end do
        k=g(i);   if(iwp1>kdiag(k))kdiag(k)=iwp1
     end if
  end do
 return
end subroutine fkdiag
subroutine invar(stress,sigm,dsbar,theta)
! forms the stress invariants in 2-d or 3-d
implicit none
  real(iwp),intent(in)::stress(:);real(iwp),intent(out)::sigm,dsbar,theta
  real(iwp)::sx,sy,sz,txy,dx,dy,dz,xj3,sine,                          &
                  s1,s2,s3,s4,s5,s6,ds1,ds2,ds3,d2,d3,sq3
  integer :: nst ; nst = ubound(stress,1)
 select case (nst)
 case(4)
  sx=stress(1); sy=stress(2); txy=stress(3); sz=stress(4)
  sigm=(sx+sy+sz)/3._iwp
  dsbar=sqrt((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+6._iwp*txy**2)/sqrt(2._iwp)
  if(dsbar<1.e-10_iwp) then
     theta=.0_iwp
  else
     dx=(2._iwp*sx-sy-sz)/3._iwp; dy=(2._iwp*sy-sz-sx)/3._iwp
     dz=(2._iwp*sz-sx-sy)/3._iwp; xj3=dx*dy*dz-dz*txy**2
     sine=-13.5_iwp*xj3/dsbar**3
     if(sine>1._iwp) sine=1._iwp
     if(sine<-1._iwp) sine=-1._iwp
     theta=asin(sine)/3._iwp
  end if
 case(6)
  sq3=sqrt(3._iwp);  s1=stress(1)  ;  s2=stress(2)
  s3=stress(3) ;  s4=stress(4);  s5=stress(5);  s6=stress(6)
  sigm=(s1+s2+s3)/3._iwp
  d2=((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/6._iwp+s4*s4+s5*s5+s6*s6
  ds1=s1-sigm ;  ds2=s2-sigm  ;  ds3=s3-sigm
  d3=ds1*ds2*ds3-ds1*s5*s5-ds2*s6*s6-ds3*s4*s4+2._iwp*s4*s5*s6
  dsbar=sq3*sqrt(d2)
  if(dsbar==0._iwp)then
      theta=0._iwp
    else
      sine=-3._iwp*sq3*d3/(2._iwp*sqrt(d2)**3)
      if(sine>1._iwp)sine=1._iwp ;  if(sine<-1._iwp)sine=-1._iwp
      theta=asin(sine)/3._iwp
  end if
 case default
  print*,"wrong size for nst in invar"
 end select
 return
end subroutine invar
subroutine formm(stress,m1,m2,m3)
! forms the derivatives of the invariants with respect to stress 2- or 3-d
 implicit none
  real(iwp),intent(in)::stress(:)
  real(iwp),intent(out)::m1(:,:),m2(:,:),m3(:,:)
  real(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm ; integer::nst , i , j
  nst=ubound(stress,1)
  select case (nst)
  case(4)
  sx=stress(1); sy=stress(2); txy=stress(3); sz=stress(4)
  dx=(2.*sx-sy-sz)/3.; dy=(2.*sy-sz-sx)/3.; dz=(2.*sz-sx-sy)/3.
  sigm=(sx+sy+sz)/3.
  m1=.0; m2=.0; m3=.0
  m1(1,1:2)=1.; m1(2,1:2)=1.; m1(4,1:2)=1.
  m1(1,4)=1.; m1(4,4)=1.; m1(2,4)=1.
  m1=m1/9./sigm
  m2(1,1)=.666666666666666; m2(2,2)=.666666666666666; m2(4,4)=.666666666666666
  m2(2,4)=-.333333333333333;m2(4,2)=-.333333333333333;m2(1,2)=-.333333333333333
  m2(2,1)=-.333333333333333;m2(1,4)=-.333333333333333;m2(4,1)=-.333333333333333
  m2(3,3)=2.; m3(3,3)=-dz
  m3(1:2,3)=txy/3.; m3(3,1:2)=txy/3.; m3(3,4)=-2.*txy/3.; m3(4,3)=-2.*txy/3.
  m3(1,1)=dx/3.; m3(2,4)=dx/3.; m3(4,2)=dx/3.
  m3(2,2)=dy/3.; m3(1,4)=dy/3.; m3(4,1)=dy/3.
  m3(4,4)=dz/3.; m3(1,2)=dz/3.; m3(2,1)=dz/3.
 case(6)
  sx=stress(1); sy=stress(2)    ;   sz=stress(3)
  txy=stress(4)  ;   tyz=stress(5) ;   tzx=stress(6)
  sigm=(sx+sy+sz)/3.
  dx=sx-sigm  ;   dy=sy-sigm ;  dz=sz-sigm
  m1 = .0; m2 = .0; m1(1:3,1:3) = 1./(3.*sigm)
  do  i=1,3 ; m2(i,i)=2. ;  m2(i+3,i+3)=6. ; end do
  m2(1,2)=-1.; m2(1,3)=-1. ; m2(2,3)=-1.; m3(1,1)=dx
  m3(1,2)=dz ; m3(1,3)=dy ; m3(1,4)=txy  ;  m3(1,5)=-2.*tyz
  m3(1,6)=tzx ; m3(2,2)=dy ; m3(2,3)=dx ; m3(2,4)=txy
  m3(2,5)=tyz ; m3(2,6)=-2.*tzx ;  m3(3,3)=dz
  m3(3,4)=-2.*txy; m3(3,5)=tyz ;  m3(3,6)=tzx
  m3(4,4)=-3.*dz ;  m3(4,5)=3.*tzx;  m3(4,6)=3.*tyz
  m3(5,5)=-3.*dx;  m3(5,6)=3.*txy ;  m3(6,6)=-3.*dy
  do  i=1,6 ;  do  j=i,6
      m1(i,j)=m1(i,j)/3.;  m1(j,i)=m1(i,j) ;  m2(i,j)=m2(i,j)/3.
      m2(j,i)=m2(i,j)   ;  m3(i,j)=m3(i,j)/3. ; m3(j,i)=m3(i,j)
  end do; end do
 case default
  print*,"wrong size for nst in formm"
 end select
 return   
 end subroutine formm              
  subroutine fmkdke(km,kp,c,ke,kd,theta)
  ! builds up 'coupled' stiffnesses ke and kd from 'elastic'
  ! stiffness km , fluid stiffness kp and coupling matrix c
  implicit none
  real(iwp),intent(in)::km(:,:),kp(:,:),c(:,:),theta
  real(iwp),intent(out)::ke(:,:),kd(:,:)
  integer::idof; idof=size(km,1)
  ke(1:idof,1:idof)=theta*km; ke(1:idof,idof+1:)=theta*c
  ke(idof+1:,1:idof)=theta*transpose(c); ke(idof+1:,idof+1:)=-theta**2*kp
  kd(1:idof,1:idof)=(theta-1.)*km; kd(1:idof,idof+1:)=(theta-1.)*c
  kd(idof+1:,1:idof)=ke(idof+1:,1:idof)
  kd(idof+1:,idof+1:)=theta*(1.-theta)*kp
  end subroutine fmkdke
subroutine banmul(kb,loads,ans)
implicit none
real(iwp),intent(in)::kb(:,:),loads(0:); real(iwp),intent(out):: ans(0:)
integer :: neq,nband,i,j ; real(iwp) :: x
neq = ubound(kb,1) ; nband = ubound(kb,2) -1
  do i = 1 , neq
     x = .0
     do j = nband + 1 , 1 , -1
        if(i+j>nband+1) x = x +kb(i,j)*loads(i+j-nband-1)
     end do
     do j = nband , 1 , -1
        if(i-j<neq - nband) x = x + kb(i-j+nband+1,j)*loads(i-j+nband+1)
     end do
     ans(i) = x
  end do
 return
end subroutine banmul
subroutine chobk1(kb,loads)
implicit none
real(iwp),intent(in)::kb(:,:); real(iwp),intent(in out)::loads(0:)
integer :: iw,n,i,j,k   ;   real(iwp) :: x
n = ubound(kb,1); iw = ubound(kb,2) - 1
  loads(1) = loads(1) / kb(1,iw+1)
    do i = 2 , n
       x = .0; k = 1
       if(i <= iw + 1 ) k = iw - i + 2   ! zero sized arrays ?
       do j = k , iw ; x=x+kb(i,j)*loads(i+j-iw-1) ; end do
       loads(i) = ( loads(i) -x)/kb(i,iw+1)
    end do
  return
end subroutine chobk1
subroutine chobk2(kb,loads)
implicit none
real(iwp),intent(in)::kb(:,:); real(iwp),intent(in out)::loads(0:)
integer :: iw,n,i,j,l,m  ;  real(iwp) :: x
 n = ubound(kb,1); iw = ubound(kb,2) - 1
 loads(n) = loads(n)/kb(n,iw+1)
 do i = n-1,1,-1
    x = .0; l = i + iw
    if(i > n-iw) l = n ; m = i + 1
    do j = m , l; x=x+kb(j,iw+i-j+1) * loads(j) ; end do
    loads(i) = (loads(i) - x)/ kb(i,iw+1)
 end do
return
end subroutine chobk2     
 subroutine fmbeam(der2,fun,points,i,ell)
!
! this subroutine forms the beam shape functions
! and their 2nd derivatives in local coordinates
!
 implicit none
 real(iwp),intent(in)::points(:,:),ell
 integer,intent(in)::i
 real(iwp),intent(out)::der2(:),fun(:)
 real(iwp)::xi,xi2,xi3
 xi=points(i,1); xi2=xi*xi; xi3=xi2*xi
 fun(1)=.25*(xi3-3.*xi+2.); fun(2)=.125*ell*(xi3-xi2-xi+1.)
 fun(3)=.25*(-xi3+3.*xi+2.); fun(4)=.125*ell*(xi3+xi2-xi-1.)
 der2(1)=1.5*xi; der2(2)=.25*ell*(3.*xi-1.)
 der2(3)=-1.5*xi; der2(4)=.25*ell*(3.*xi+1.)
 end subroutine fmbeam                                                       
 subroutine beam_km(km,ei,ell)
!
!      this subroutine forms the stiffness matrix of a
!      line beam element(bending only)
!
 implicit none
   real(iwp),intent(in)::ei,ell; real,intent(out):: km(:,:)
   km(1,1)=12.*ei/(ell*ell*ell) ;km(3,3)=km(1,1)
   km(1,2)=6.*ei/(ell*ell) ;  km(2,1)=km(1,2) ; km(1,4)=km(1,2)
   km(4,1)=km(1,4) ;  km(1,3)=-km(1,1) ; km(3,1)=km(1,3) ; km(3,4)=-km(1,2)
   km(4,3)=km(3,4) ; km(2,3)=km(3,4) ;  km(3,2)=km(2,3); km(2,2)=4.*ei/ell
   km(4,4)=km(2,2) ;km(2,4)=2.*ei/ell ; km(4,2)=km(2,4)
  return
 end subroutine beam_km                                                      
 subroutine beam_mm(mm,fs,ell)
 implicit none
 real(iwp),intent(in)::fs,ell
 real(iwp),intent(out)::mm(:,:)
 real(iwp)::fac
!
!     this subroutine forms the mass matrix of a beam element
!
 fac=(fs*ell)/420. 
 mm(1,1)=156.*fac; mm(3,3)=mm(1,1); mm(1,2)=22.*ell*fac; mm(2,1)=mm(1,2)
 mm(3,4)=-mm(1,2); mm(4,3)=mm(3,4); mm(1,3)=54.*fac; mm(3,1)=mm(1,3)
 mm(1,4)=-13.*ell*fac; mm(4,1)=mm(1,4); mm(2,3)=-mm(1,4)
 mm(3,2)=mm(2,3); mm(2,2)=4.*(ell**2)*fac; mm(4,4)=mm(2,2)
 mm(2,4)=-3.*(ell**2)*fac; mm(4,2)=mm(2,4)
 return
 end subroutine beam_mm  
subroutine beam_kp(kp,coord,pax)
!
!      this subroutine forms the terms of the beam stiffness
!      matrix due to axial loading
!
implicit none
real(iwp),intent(in)::coord(:,:),pax; real(iwp),intent(out)::kp(:,:)
real(iwp)::x1,x2,y1,y2,ell,c,s,pre(6,4),k_1d(4,4)
integer::ndim
 ndim=ubound(coord,2)
 select case(ndim)
 case(1)
 x1=coord(1,1); x2=coord(2,1)
 ell=x2-x1
 kp(1,1)=1.2/ell; kp(1,2)=0.1
 kp(2,1)=0.1; kp(1,3)=-1.2/ell
 kp(3,1)=-1.2/ell; kp(1,4)=0.1
 kp(4,1)=0.1; kp(2,2)=2.0*ell/15.0
 kp(2,3)=-0.1; kp(3,2)=-0.1
 kp(2,4)=-ell/30.0; kp(4,2)=-ell/30.0
 kp(3,3)=1.2/ell; kp(3,4)=-0.1
 kp(4,3)=-0.1; kp(4,4)=2.0*ell/15.0
 kp=kp*pax
 case(2)
 x1=coord(1,1); y1=coord(1,2) 
 x2=coord(2,1); y2=coord(2,2)
 ell=sqrt((y2-y1)**2+(x2-x1)**2)
 c=(x2-x1)/ell; s=(y2-y1)/ell
 pre=0.0
 pre(1,1)=-s; pre(2,1)=c; pre(3,2)=1.0; pre(4,3)=-s; pre(5,3)=c; pre(6,4)=1.0
 k_1d(1,1)=1.2/ell;   k_1d(1,2)=0.1
 k_1d(2,1)=0.1;       k_1d(1,3)=-1.2/ell
 k_1d(3,1)=-1.2/ell;  k_1d(1,4)=0.1
 k_1d(4,1)=0.1;       k_1d(2,2)=2.0*ell/15.0
 k_1d(2,3)=-0.1;      k_1d(3,2)=-0.1
 k_1d(2,4)=-ell/30.0; k_1d(4,2)=-ell/30.0
 k_1d(3,3)=1.2/ell;   k_1d(3,4)=-0.1
 k_1d(4,3)=-0.1;      k_1d(4,4)=2.0*ell/15.0
 kp=matmul(matmul(pre,k_1d),transpose(pre))
 kp=kp*pax
 end select
 return
end subroutine beam_kp   

 subroutine interp(k,dt,rt,rl,al,ntp)
!
!     this subroutine forms the load/time functions by interpolation
!     if dt is not an exact multiple it stops one short
!
 real(iwp),intent(in)::dt,rt(:),rl(:)
 integer,intent(in)::k,ntp; real(iwp),intent(out)::al(:,:)
 integer::np,i,j; real(iwp)::t,val 
 np=size(rt); al(1,k)=rl(1); t=rt(1)
 do j=2,ntp
   t=t+dt
   do i=2,np
     if(t.le.rt(i))then
       val=rl(i-1)+((t-rt(i-1))/(rt(i)-rt(i-1)))*(rl(i)-rl(i-1))
       exit
     end if
   end do
   al(j,k)=val
 end do
 return
 end subroutine interp
                                                                             

 subroutine find_g(num,g,rest)
 ! Finds g from node numbers and restraints "rest"
 ! Only for uncoupled problems
 implicit none
 integer,intent(in)::num(:),rest(:,:); integer,intent(out)::g(:)
 integer::inc,i,j,k,l,count,nr,nod,nodof
 nr=ubound(rest,1); nod=ubound(num,1);nodof=ubound(rest,2)-1
 g=1;inc=0
 do i=1,nod
  count=0;l=num(i)
  do j=1,nr
     if(l>rest(j,1)) then
        do k=2,nodof+1; if(rest(j,k)==0)count=count+1;end do
     end if
  end do
  do k=2,nodof+1
   inc=inc+1
   do j=1,nr
    if(l==rest(j,1).and.rest(j,k)==0) then
      g(inc)=0; count=count+1
    end if
   end do
   if(g(inc)/=0) g(inc) = l*nodof-count-(nodof+1-k)
  end do
 end do
 end subroutine find_g

subroutine find_g2(num,g,rest)
! second effort at find_g
 implicit none ; integer,intent(in)::num(:),rest(:,:)
 integer,intent(out)::g(:); integer::i,j,l,s1,s2,nr,nodof,nod
 nr = ubound(rest,1); nodof = ubound(rest,2)-1; nod = ubound(num,1)
   DO i=1,nod ; l= num(i)
      DO j=1,nr
       IF(l/=rest(j,1)) THEN
        IF(l>rest(j,1)) THEN; s1=rest(j,1); s2=MAXVAL(rest(j,2:)); END IF
        END IF
      END DO
      g(3*i-2) = s2+3*(l-s1)-2; g(3*i-1) = s2+3*(l-s1)-1; g(3*i) = s2+3*(l-s1)
   END DO
   DO i=1,nod;   l = num(i)
     DO j = 1 , nr
        IF(l==rest(j,1)) THEN
            g(3*i-2)=rest(j,2); g(3*i-1)=rest(j,3); g(3*i)=rest(j,4) ; EXIT
        END IF
     END DO
   END DO
end subroutine find_g2
subroutine rearrange(rest)
 implicit none; integer,intent(inout)::rest(:,:)
 integer::nr,nodof,i,k,m; nr=ubound(rest,1); nodof = ubound(rest,2)-1
  m=0
   DO i=1,nr; DO k=2,nodof+1
      IF(rest(i,k)/=0) THEN; m=m+1; rest(i,k)=m; END IF
   END DO; END DO
   DO i=1,nr; k=rest(i,1)
    DO m=2,nodof+1; IF(rest(i,m)/=0)rest(i,m) = rest(i,m) + nodof*(k-i); END DO
   END DO
end subroutine rearrange

subroutine find_g3(num,g,rest)
! 3rd try - binary search
integer,intent(in)::rest(:,:),num(:); integer,intent(out)::g(:)
integer:: i,j,l,nod,nodof,s1,s2,s3,only,first,last,half
logical::found       ; nod = ubound(num,1) ; nodof = ubound(rest,2) - 1
do i=1,nod ; l = num(i)
 first = 1 ;  last = ubound(rest,1)
  do
   if(first==last) exit ! silly array
   half = (first + last)/2
   if(l<=rest(half,1)) then
     last = half  ! discard second half
   else
     first = half + 1  ! discard first half
   end if
  end do
 only = first; found = (l==rest(only,1))
 IF(found) THEN
  DO j = 1 , nodof
   g(nodof*i-nodof+j) = rest(only,j+1)  
  END DO
 ELSE
  s1 = only -1 ;   s2 = maxval(rest(s1,2:))    ; s3 = l - rest(s1,1)
  DO j=1,nodof
   g(nodof*i-nodof+j) = s2+nodof*s3-nodof+j
  END DO
 END IF
end do
end subroutine find_g3

subroutine rearrange_2(rest)
 implicit none; integer,intent(inout)::rest(:,:)
 integer::nr,i,m; nr=ubound(rest,1)
   DO i=1,nr ; rest(i,2) = rest(i,1) - i; END DO
   DO i=1,nr;  IF(rest(i,2)/=0) THEN ; m=i-1; EXIT ; END IF; END DO
   DO i=m,nr;  rest(i,2) = rest(i,2) + 1; END DO
end subroutine rearrange_2

subroutine find_g4(num,g,rest)
! 4th try - binary search  for one degree of freedom per node
integer,intent(in)::rest(:,:),num(:); integer,intent(out)::g(:)
integer:: i,l,nod,only,first,last,half ,s1,s2,s3
logical::found       ; nod = ubound(num,1) ; nodof = ubound(rest,2) - 1
do i=1,nod ; l = num(i)
 first = 1 ;  last = ubound(rest,1)
  do
   if(first==last) exit ! silly array
   half = (first + last)/2
   if(l<=rest(half,1)) then
     last = half  ! discard second half
   else
     first = half + 1  ! discard first half
   end if
  end do
 only = first; found = (l==rest(only,1))
 IF(found) THEN ;  g(i) = 0
 ELSE ; s1 = only -1; s2 = rest(s1,2); s3 = l - rest(s1,1)
  g(i) = s2 + s3 - 1
 END IF
end do
end subroutine find_g4

subroutine bandred(a,d,e)
!
! This subroutine transforms a real symmetric band matrix a,
! of order n and band width iw, to tridiagonal form by an appropriate
! sequence of Jacobi rotations. During the transformation the
! property of the band matrix is maintained. The method yields
! a tridiagonal matrix, the diagonal elements of which are in
! d(n) and off-diagonal elements in e(n).
!
implicit none  
real(iwp),intent(in out)::a(:,:)
real(iwp),intent(out)::d(0:),e(0:)
integer::  iw, n2, n, k, maxr, irr, ir, kr, j, jm, iugl, j2,    &
           l, jl, maxl, i
real(iwp) :: g, b, s, c, c2, s2, cs, u, u1
n=ubound(a,1)
iw = ubound(a,2)-1
n2 = n - 2
if (n2>=1) then
  loop_1: do k=1,n2
    maxr = iw
    if (n-k<iw) maxr = n - k
    loop_2: do irr=2,maxr
      ir = 2 + maxr - irr
      kr = k + ir
      loop_3: do j=kr,n,iw
        if (j/=kr) then
          if (abs(g) < tiny(g)) exit loop_3
          jm = j - iw
          b = -a(jm-1,iw+1)/g
          iugl = j - iw
        else
          if (abs(a(k,ir+1))<tiny(a(k,ir+1))) exit loop_3
          b = -a(k,ir)/a(k,ir+1)
          iugl = k
        end if
        s = 1./sqrt(1.+b*b)
        c = b*s
        c2 = c*c
        s2 = s*s
        cs = c*s
        u = c2*a(j-1,1) - 2.*cs*a(j-1,2) + s2*a(j,1)
        u1 = s2*a(j-1,1) + 2.*cs*a(j-1,2) + c2*a(j,1)
        a(j-1,2) = cs*(a(j-1,1)-a(j,1)) + (c2-s2)*a(j-1,2)
        a(j-1,1) = u
        a(j,1) = u1
        j2 = j - 2
        do  l=iugl,j2
          jl = j - l
          u = c*a(l,jl) - s*a(l,jl+1)
          a(l,jl+1) = s*a(l,jl) + c*a(l,jl+1)
          a(l,jl) = u
        end do
        jm = j - iw
        if (j/=kr) a(jm-1,iw+1) = c*a(jm-1,iw+1) - s*g
        maxl = iw - 1
        if (n-j<iw-1) maxl = n - j
        if (maxl>0) then
          do  l=1,maxl
            u = c*a(j-1,l+2) - s*a(j,l+1)
            a(j,l+1) = s*a(j-1,l+2) + c*a(j,l+1)
            a(j-1,l+2) = u
          end do
        end if
        if (j+iw<=n) then
          g = -s*a(j,iw+1)
          a(j,iw+1) = c*a(j,iw+1)
        end if
      end do loop_3
    end do loop_2
  end do loop_1
end if
e(1) = 0.
d(1:n) = a(1:n,1)
if (2<=n) then
  do i=2,n
    e(i) = a(i-1,2)
  end do
end if
return
end subroutine bandred

subroutine bisect(d,e,acheps,ifail)
!
! This subroutine finds the eigenvalues of a tridiagonal matrix,
! given with its diagonal elements in the array d(n) and
! its subdiagonal elements in the last n - 1 stores of the
! array e(n), using ql transformations. The eigenvalues are
! overwritten on the diagonal elements in the array d in
! ascending order. The subroutine will fail if any one
! eigenvalue takes more than 30 iterations.
!
implicit none
real(iwp),intent(in)::acheps
real(iwp),intent(in out)::d(0:),e(0:)
integer,intent(out)::ifail  
integer :: m, n, i, l, j, i1, m1, ii, aux
real(iwp) :: b, f, h, g, p, r, c, s
n = ubound(d,1)
if (n/=1) then
  do  i=2,n
    e(i-1) = e(i)
  end do
end if
e(n) = 0.
b = 0.
f = 0.
loop_1: do l=1,n
  j = 0
  h = acheps*(abs(d(l))+abs(e(l)))
  if (b<h) b = h
  ! look for small sub diagonal element
  do  m=l,n
    if (abs(e(m))<=b) exit
  end do
  if (m/=l) then
    loop_2: do
      if (j==30) then
        ifail=1
        exit loop_1
      end if
      j = j + 1
      ! form shift
      g = d(l)
      h = d(l+1) - g
      if (abs(h) < abs(e(l))) then
        p = h*0.5/e(l)
        r = sqrt(p*p+1.)
        h = p + r
        if (p<0.) h = p - r
        d(l) = e(l)/h
      else
        p = 2.*e(l)/h
        r = sqrt(p*p+1.)
        d(l) = e(l)*p/(1.+r)
      end if
      h = g - d(l)
      i1 = l + 1
      if (i1<=n) then
        do i=i1,n
          d(i) = d(i) - h
        end do
      end if
      f = f + h
      ! ql transformation
      p = d(m)
      c = 1.
      s = 0.
      m1 = m - 1
      loop_4: do ii=l,m1
        i = m1 - ii + l
        g = c*e(i)
        h = c*p
        if (abs(p)>=abs(e(i))) then
          c = e(i)/p
          r = sqrt(c*c+1.)
          e(i+1) = s*p*r
          s = c/r
          c = 1./r
        else
          c = p/e(i)
          r = sqrt(c*c+1.)
          e(i+1) = s*e(i)*r
          s = 1./r
          c = c/r
        end if
        p = c*d(i) - s*g
        d(i+1) = h + s*(c*g+s*d(i))
      end do loop_4
      e(l)= s*p
      d(l)= c*p
      if (abs(e(l))<=b) exit loop_2
    end do loop_2
  end if
  p = d(l) + f
  ! order eigenvalue
  aux=0
  if (l/=1) then
    loop_3: do  ii=2,l
      i = l - ii + 2
      if (p>=d(i-1)) then
        aux=1
        exit loop_3
      end if
      d(i) = d(i-1)
    end do loop_3
  end if
  if(aux==0)then
    i = 1
  end if
  d(i) = p
  ifail=0
end do loop_1
return
end subroutine bisect
end module new_library
