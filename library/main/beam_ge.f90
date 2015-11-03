SUBROUTINE beam_ge(ge,ell)
!
! This subroutine forms the beam geometric matrix for stability analysis.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::ell
 REAL(iwp),INTENT(OUT)::ge(:,:)
 REAL(iwp)::pt1=0.1_iwp,opt2=1.2_iwp,two=2.0_iwp,d15=15.0_iwp,d30=30.0_iwp
 ge(1,1)=opt2/ell
 ge(1,2)=pt1
 ge(2,1)=pt1
 ge(1,3)=-opt2/ell
 ge(3,1)=-opt2/ell
 ge(1,4)=pt1
 ge(4,1)=pt1
 ge(2,2)=two*ell/d15
 ge(2,3)=-pt1
 ge(3,2)=-pt1
 ge(2,4)=-ell/d30
 ge(4,2)=-ell/d30
 ge(3,3)=opt2/ell
 ge(3,4)=-pt1
 ge(4,3)=-pt1
 ge(4,4)=two*ell/d15
RETURN
END SUBROUTINE beam_ge
