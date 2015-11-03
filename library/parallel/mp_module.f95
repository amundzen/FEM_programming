MODULE mp_module
!
! Variables needed for message passing version of code.
!
  USE precision
  IMPLICIT NONE
!
!     Include file(s):
!


  INCLUDE "mpif.h"

!
!     Variables for packing data to be broadcast.
!

  integer :: bufsizer, position, bufsize, bufdecl, recbufsize

  integer, PARAMETER :: ilength = 8, rlength = 8


!

  INTEGER, ALLOCATABLE :: tempbuf(:)

!#ifdef sp2
!      CHARACTER(LEN=1), ALLOCATABLE :: tempbuf(:)
!#endif


  integer :: mytid, ier
  integer :: ielpe, ic_local
!
!     Variables used only for identification of columns associated with each PE
!

  integer :: ggstart, ggend, vecstart, vecend
  REAL(iwp) :: sum_load_temp, dotrdold_temp
  integer, ALLOCATABLE :: nf_temp(:,:)
  REAL(iwp), ALLOCATABLE :: totd_temp2(:)
  REAL(iwp), ALLOCATABLE :: g_coord_temp(:,:)

END MODULE mp_module
