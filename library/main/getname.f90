SUBROUTINE getname(argv,nlen)
!
! This subroutine the base name of data file.
!
 IMPLICIT NONE
 INTEGER::narg
 INTEGER,INTENT(OUT)::nlen
 INTEGER::lnblnk,iargc
 CHARACTER(*),INTENT(OUT)::argv
 LOGICAL found
 narg=iargc()
 IF(narg.lt.1)THEN
   WRITE(*,*)'Please enter the base name of data file: '
   READ(*,*) argv
  ELSE
   CALL getarg(1,argv)
 ENDIF
 nlen=lnblnk(argv)
 INQUIRE(file=argv(1:nlen)//'.dat',exist=found)
 IF(.not.found)THEN
  WRITE(*,*)'Data file not found: ',argv(1:nlen)//'.dat'
  WRITE(*,*)'Please create or check spelling.'
  STOP
 ENDIF
RETURN
END SUBROUTINE getname
