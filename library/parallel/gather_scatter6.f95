MODULE gather_scatter6
!
! Version 1c, M.A. Pettipher, Manchester Computing: 9th December 1996
!             Combines all routines associated with gather/scatter code
!             into one module. 
!             Removes references to tempbuf, bufsize, bufdecl.
!             Updated from Version 1a 24th July 2001 - added scatter2
!             Updated from Version 1a 3rd September 2001 - added scatter3
!
  USE precision
  USE global_variables1
  USE utility
  IMPLICIT NONE
!
! Variables restricted to gather_scatter module:
!
  INTEGER, ALLOCATABLE :: ggl_pp(:,:)
  INTEGER, ALLOCATABLE, PRIVATE :: lenget(:), lenput(:), lengetsum(:)
  INTEGER, ALLOCATABLE, PRIVATE :: toget(:,:), toput(:,:)
  INTEGER, ALLOCATABLE, PRIVATE :: pesget(:), pesput(:)
  INTEGER, ALLOCATABLE, PRIVATE :: getpes(:), putpes(:)
  INTEGER, ALLOCATABLE, PRIVATE :: toget_temp(:,:), toput_temp(:,:) ! Used only in make_ggl and allocate_gather_scatter
  REAL(iwp), ALLOCATABLE, PRIVATE :: tempget(:), tempput(:), tempput1(:,:)
  REAL(iwp), ALLOCATABLE, PRIVATE :: sumget(:,:)
  INTEGER, PRIVATE :: numpesget, numpesput, numpesgetput
  INTEGER, PRIVATE :: npes_pp, len_pl_pp
!
! Variables also used in main program:
!
  INTEGER :: npes, nels_pp1, nels_pp2, neq_pp1, neq_pp2
  INTEGER :: num_nels_pp1, num_neq_pp1, threshold
  INTEGER :: ibar
  INTEGER :: iel_start, ieq_start
CONTAINS
FUNCTION calc_ielpe(iel) RESULT(ielpe)
IMPLICIT NONE
  INTEGER :: ielpe
  INTEGER, INTENT(IN) :: iel
  IF (numpe <= num_nels_pp1 .OR. num_nels_pp1 == 0) THEN
    ielpe = (numpe - 1)*nels_pp1 + iel
  ELSE
    ielpe = num_nels_pp1*nels_pp1 + (numpe - num_nels_pp1 - 1)*(nels_pp1 - 1) + iel
  ENDIF
END FUNCTION calc_ielpe
SUBROUTINE calc_npes_pp
!
! Set value to be used in initial dimensions of toget and toput
! Do not know exact number of processors needed by each processor
! until make_ggl called, but must declare these arrays (or at least
! temporary versions of them) before exact number is known. npes
! is definitely enough, but wasteful of memory, so make rough
! overestimate based on number of processors, npes.
!
IMPLICIT NONE
  SELECT CASE (npes)
    CASE (1)
      npes_pp = 1
    CASE (2)
      npes_pp = 2
    CASE (3:15)
      npes_pp = 3
    CASE DEFAULT
      npes_pp = npes/4
  END SELECT
END SUBROUTINE calc_npes_pp
SUBROUTINE calc_nels_pp
IMPLICIT NONE
!
! Calculate number of elements per processor (number of columns of
! matrix): nels_pp, nels_pp1, nels_pp2
! Also calculate first and last element numbers for each processor,
! iel_start, iel_end
!
  IF (npes == 1) THEN
    nels_pp1 = 0
    nels_pp2 = nels
    nels_pp = nels_pp2
    iel_start = 1
!   iel_end = nels
  ELSE
    nels_pp2 = nels/npes
    num_nels_pp1 = nels - nels_pp2*npes
    IF (num_nels_pp1 == 0) THEN
      nels_pp1 = nels_pp2
    ELSE
      nels_pp1 = nels_pp2 + 1
    ENDIF
    IF (numpe <= num_nels_pp1 .OR. num_nels_pp1 == 0) THEN
      nels_pp = nels_pp1
      iel_start = (numpe - 1)*nels_pp1 + 1
    ELSE
      nels_pp = nels_pp2
      iel_start = num_nels_pp1*nels_pp1 + (numpe - num_nels_pp1 - 1)*(nels_pp1 - 1) + 1
    ENDIF
!   iel_end = iel_start + nels_pp
  ENDIF
  WRITE(details,'(A,I6,A,I6)') 'PE no: ', numpe, ' nels_pp: ', nels_pp
END SUBROUTINE calc_nels_pp
SUBROUTINE calc_neq_pp
IMPLICIT NONE
!
! Calculate number of equations per processor (number of 'elements' of
! vectors): neq_pp, neq_pp1, neq_pp2
! Also calculate first and last equation numbers for each processor,
! ieq_start
!
  IF (npes == 1) THEN
    neq_pp2 = neq
    neq_pp1 = neq
    neq_pp = neq_pp1
    ieq_start = 1
  ELSE
    neq_pp2 = neq/npes
    num_neq_pp1 = neq - neq_pp2*npes
    IF (num_neq_pp1 == 0) THEN
      neq_pp1 = neq_pp2
    ELSE
      neq_pp1 = neq_pp2 + 1
    ENDIF
    IF (numpe <= num_neq_pp1 .OR. num_neq_pp1 == 0) THEN
      neq_pp = neq_pp1
      ieq_start = (numpe - 1)*neq_pp1 + 1
    ELSE
      neq_pp = neq_pp2
      ieq_start = num_neq_pp1*neq_pp1 + (numpe - num_neq_pp1 - 1)*(neq_pp1 - 1) + 1
    ENDIF
  ENDIF
  WRITE(details,'(A,I6,A,I6)') 'PE no: ', numpe, ' neq_pp: ', neq_pp
END SUBROUTINE calc_neq_pp
SUBROUTINE allocate_gather_scatter
!
! Allocates arrays needed in gather_scatter module
!
  ALLOCATE ( ggl_pp(ntot,nels_pp) )
  ALLOCATE ( lenget(npes), lenput(npes), lengetsum(0:npes) )
  ALLOCATE ( pesget(npes), pesput(npes) )
  ALLOCATE ( getpes(npes), putpes(npes) )
!
! Not obvious what are reasonable values for array dimensions. With 1000 element problem
! and 32 PEs, neq_pp1 = 413, but PE 31 requires 832 remote accesses, so array size must be
! more than 2*neq_pp1 - so set to 3*neq_pp1.
!
  ALLOCATE ( tempget(3*neq_pp1), tempput(3*neq_pp1), tempput1(3*neq_pp1,npes))
!
! Allocate temporary arrays for toget and toput, until actual size
! is known - reallocate later.
!
  ALLOCATE ( toget_temp(neq_pp1,npes_pp), toput_temp(neq_pp1,npes_pp) )
!
! Allocate sumget with npes_pp as well
!
  ALLOCATE ( sumget(neq_pp1,npes_pp) )
END SUBROUTINE allocate_gather_scatter
SUBROUTINE deallocate_gather_scatter
!
! Deallocates arrays needed in gather_scatter module
!
  DEALLOCATE ( ggl_pp )
  DEALLOCATE ( lenget, lenput, lengetsum )
! DEALLOCATE ( toget, toput )
  DEALLOCATE ( pesget, pesput )
  DEALLOCATE ( getpes, putpes )
  DEALLOCATE ( tempget, tempput, tempput1 )
END SUBROUTINE deallocate_gather_scatter
!
SUBROUTINE gather(p_pp, pmul_pp)
!
!     Performs the gather operation: pmul = p(g)
!
  IMPLICIT NONE
!
!     Include file(s):
!
  INCLUDE "mpif.h"
!
!     Declare arguments:
!
  REAL(iwp), INTENT(IN) :: p_pp(neq_pp)
  REAL(iwp), INTENT(OUT) :: pmul_pp(ntot,nels_pp)
!     Local Variables:
      INTEGER i, ii, j, k, ier, nstart, pe_number
  INTEGER vrequest(numpesput)
  INTEGER vstatus(MPI_STATUS_SIZE, numpesgetput)
  REAL(iwp), DIMENSION(neq_pp1,numpesput) :: vtempput
  LOGICAL lflag
  INTEGER recbufsize
  INTEGER status(MPI_STATUS_SIZE)
  REAL(iwp), DIMENSION(0:len_pl_pp) :: pl_pp
!
! Barrier to synchronise before gathering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 20
       CALL my_barrier(numpe,ibar,details,'First barrier in gather')
! *****
! *****   (B1) - Calculate local p vector, pl_pp
! *****   
  pl_pp = 0.0
  DO ii = 1, numpesput
    i = pesput(ii)
      DO j = 1, lenput(i)
        vtempput(j,ii) = p_pp(toput(j,ii))
      END DO
      CALL MPI_ISEND (vtempput(1,ii),lenput(i),MPI_REAL8,i-1,numpe,           &
        MPI_COMM_WORLD,vrequest(ii),ier)
      IF (ier .NE. 0) CALL mperror('Error in (B1) send',ier)
!---- write number of processes to gather from
!      WRITE(11,*) 'Processor ',numpe,' gathers from ',numpesput,' processors'
  END DO
! 
! Each PE must obtain required values of p
!
! 
! For data already available, just make a local copy:
!
  nstart = lengetsum(numpe-1)
  DO j = 1, lenget(numpe)
    pl_pp(nstart + j) = p_pp(toget(j,getpes(numpe)))
  END DO
! 
! For non local data must receive it:
!
  DO i = 1, numpesget
    CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,                       &
      MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B1) probe',ier)
    CALL MPI_GET_COUNT (vstatus(1,i),MPI_REAL8,recbufsize,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B1) get_count',ier)
    pe_number = vstatus(MPI_tag,i)

    nstart = lengetsum(pe_number-1)
    CALL MPI_RECV (pl_pp(nstart + 1),lenget(pe_number),MPI_REAL8,                    &
!     MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,i),ier)
      vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B1) receive',ier)
  END DO
! CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
! IF (ier .NE. 0) CALL mperror('Error in (B1) barrier after receives',ier)          
! *****
! *****   (B2) - Generate local pmul, pmul_pp
! *****
  DO j = 1,nels_pp
    pmul_pp(:,j) = pl_pp(ggl_pp(:,j))
  END DO
!
! Make sure all ISENDs have completed and free up internal book-keeping of requests
!
  CALL MPI_WAITALL(numpesput,vrequest,vstatus,ier)
  IF (ier /= MPI_SUCCESS) CALL mperror ('Error in MPI_WAITALL', ier)
!
! Barrier to synchronise after gathering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 21
       CALL my_barrier(numpe,ibar,details,'Last barrier in gather')
END SUBROUTINE gather
SUBROUTINE gather1(p_pp, pmul_pp)
!
!     Performs the gather operation: pmul = p(g)
!
  IMPLICIT NONE
!
!     Include file(s):
!
  INCLUDE "mpif.h"
!
!     Declare arguments:
!  
  REAL(iwp), INTENT(IN) :: p_pp(neq_pp)
  REAL(iwp), INTENT(OUT) :: pmul_pp(ntot,nels_pp)
!     Local Variables:
      INTEGER i, ii, j, k, ier, nstart, pe_number
  INTEGER vrequest(numpesput)
  INTEGER vstatus(MPI_STATUS_SIZE, numpesgetput)
  REAL(iwp), DIMENSION(neq_pp1,numpesput) :: vtempput
  LOGICAL lflag
      INTEGER recbufsize
      INTEGER status(MPI_STATUS_SIZE)
  REAL(iwp), DIMENSION(0:len_pl_pp) :: pl_pp
!
! Barrier to synchronise before gathering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
!      ibar = 20
!      CALL my_barrier(numpe,ibar,details,'First barrier in gather')
! *****
! *****   (B1) - Calculate local p vector, pl_pp
! *****   
  pl_pp = 0.0
  DO ii = 1, numpesput
    i = pesput(ii)
      DO j = 1, lenput(i)
        vtempput(j,ii) = p_pp(toput(j,ii))
      END DO
      CALL MPI_ISEND (vtempput(1,ii),lenput(i),MPI_REAL8,i-1,numpe,           &
        MPI_COMM_WORLD,vrequest(ii),ier)
      IF (ier .NE. 0) CALL mperror('Error in (B1) send',ier)
!---- write number of processes to gather from
!      WRITE(11,*) 'Processor ',numpe,' gathers from ',numpesput,' processors'
  END DO
! 
! Each PE must obtain required values of p
!
! 
! For data already available, just make a local copy:
!
  nstart = lengetsum(numpe-1)
  DO j = 1, lenget(numpe)
    pl_pp(nstart + j) = p_pp(toget(j,getpes(numpe)))
  END DO
! 
! For non local data must receive it:
!
  DO i = 1, numpesget
    CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,                       &
      MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B1) probe',ier)
    CALL MPI_GET_COUNT (vstatus(1,i),MPI_REAL8,recbufsize,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B1) get_count',ier)
    pe_number = vstatus(MPI_tag,i)
    nstart = lengetsum(pe_number-1)
    CALL MPI_RECV (pl_pp(nstart + 1),lenget(pe_number),MPI_REAL8,                    &
!     MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,i),ier)
      vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B1) receive',ier)
  END DO
! CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
! IF (ier .NE. 0) CALL mperror('Error in (B1) barrier after receives',ier)            
! *****
! *****   (B2) - Generate local pmul, pmul_pp
! *****
  DO j = 1,nels_pp
    pmul_pp(:,j) = pl_pp(ggl_pp(:,j))
  END DO
!
! Make sure all ISENDs have completed and free up internal book-keeping of requests
!
  CALL MPI_WAITALL(numpesput,vrequest,vstatus,ier)
  IF (ier /= MPI_SUCCESS) CALL mperror ('Error in MPI_WAITALL', ier)
!
! Barrier to synchronise after gathering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 21
       CALL my_barrier(numpe,ibar,details,'Last barrier in gather')
END SUBROUTINE gather1
SUBROUTINE scatter(u_pp, utemp_pp)
!
!     Performs the scatter operation: u(g) = u(g) + utemp
!
  IMPLICIT NONE
!
! Include file(s):
!
  INCLUDE "mpif.h" 
!
!     Arguments:
!
  REAL(iwp), INTENT(INOUT) :: u_pp(neq_pp)
  REAL(iwp) :: utemp_pp(ntot,nels_pp)            ! not needed as argument(use automatic array)
  INTEGER vrequest(numpesget)
  INTEGER vstatus(MPI_STATUS_SIZE, numpesgetput)
!
!     Local Variables:
!
  INTEGER i, j, k, ii, ier, nstart, pe_number
  INTEGER ibar
  INTEGER recbufsize
  INTEGER status(MPI_STATUS_SIZE)
  INTEGER temp_status(MPI_STATUS_SIZE)
  LOGICAL lflag
  REAL(iwp), DIMENSION(0:len_pl_pp) :: ul_pp
!
! Barrier to synchronise before scattering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 23
       CALL my_barrier(numpe,ibar,details,'First barrier in gather')
!
! (B4) - Generate ul_pp vector
! 
  ul_pp = 0.0_iwp
  DO j = 1,nels_pp
    ul_pp(ggl_pp(:,j)) = ul_pp(ggl_pp(:,j)) + utemp_pp(:,j)
  END DO
! DO j = 1,nels_pp
!   DO i = 1, ntot
!     ul_pp(ggl_pp(i,j)) = ul_pp(ggl_pp(i,j))           &
!       + utemp_pp(i,j)
!   END DO
! END DO
!
! (B5) - Generate u_pp vector - using ul_pp vectors from all
!        appropriate PEs
! 
!
! Local data:
!
  nstart = lengetsum(numpe-1)
  DO j = 1, lenget(numpe)
    u_pp(toget(j,getpes(numpe))) = u_pp(toget(j,getpes(numpe))) +        &
      ul_pp(nstart + j)
  END DO
! 
! Remote data:
! 
  DO ii = 1, numpesget
    i = pesget(ii)
    nstart = lengetsum(i-1)
    CALL MPI_ISEND (ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,             &
      MPI_COMM_WORLD,vrequest(ii),ier)
!   CALL MPI_SEND (ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,             &
!     MPI_COMM_WORLD,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) isend',ier)
  END DO
! *****
! *****    Now receive data
! *****
! ibar = 500
! CALL my_barrier(numpe,ibar,details,'Scatter: Barrier before receives')
! CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,     &
!   temp_status,ier)
! IF (ier .NE. 0) CALL mperror('Error in (temp) probe',ier)
  DO i = 1, numpesput
    CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,     &
      vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) probe',ier)
    CALL MPI_GET_COUNT (vstatus(1,i),MPI_REAL8,recbufsize,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) get_count',ier)
    pe_number = vstatus(MPI_TAG,i)
!#ifdef debug
    IF (recbufsize .NE. lenput(pe_number)) THEN
      WRITE(*,*) 'PE: ', numpe, ' Scatter error: recbufsize: ', recbufsize, ' /= lenput(pe_number): ', lenput(pe_number)
      WRITE(*,*) ' pe_number: ', pe_number, ' i: ', i
      WRITE(*,*) 'source: ', vstatus(MPI_SOURCE,i), ' tag: ', vstatus(MPI_TAG,i)
      CALL MPI_GET_COUNT (temp_status,MPI_REAL8,recbufsize,ier)
      IF (ier .NE. 0) CALL mperror('Error in (B5) get_count',ier)
      WRITE(*,*) 'Temp probe. recbufsize: ', recbufsize, ' source: ', temp_status(MPI_SOURCE),     &
        ' tag: ', temp_status(MPI_TAG)
      STOP
    ENDIF
!#endif
    CALL MPI_RECV (tempput(1),lenput(pe_number),MPI_REAL8,                 &
!     MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,i),ier)
      vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) receive',ier)
!
!   Note - need to use tempput to receive data, as must be added to
!   u_pp in loop below. No problem using the same temporary buffer
!   for the receives as these are blocking receives.
!
    DO j = 1, lenput(pe_number)
      u_pp(toput(j,putpes(pe_number))) = u_pp(toput(j,putpes(pe_number))) + tempput(j)
    END DO
  END DO
!
! Make sure all ISENDs have completed and free up internal book-keeping of requests
!
  CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
  IF (ier /= MPI_SUCCESS) CALL mperror ('Error in MPI_WAITALL', ier)
!
! Barrier to synchronise after scattering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 24
       CALL my_barrier(numpe,ibar,details,'Last barrier in scatter')
END SUBROUTINE scatter
SUBROUTINE scatter1(u_pp, utemp_pp)
!
!     Performs the scatter operation: u(g) = u(g) + utemp
!     Modified for setting up diag_precon: utemp is a vector
!
  IMPLICIT NONE
!
! Include file(s):
!
  INCLUDE "mpif.h" 
!
!     Arguments:
!
  REAL(iwp), INTENT(INOUT) :: u_pp(neq_pp)
  REAL(iwp) :: utemp_pp(ntot)            ! not needed as argument(use automatic array)
  INTEGER vrequest(numpesget)
  INTEGER vstatus(MPI_STATUS_SIZE, numpesget)
!
!     Local Variables:
!
  INTEGER i, j, k, ii, ier, nstart, pe_number
  INTEGER position, recbufsize
  INTEGER status(MPI_STATUS_SIZE)
  INTEGER temp_status(MPI_STATUS_SIZE)
  REAL(iwp), DIMENSION(neq_pp*4) :: ul_pp
  REAL(iwp) :: tempputs(neq_pp1,npes)
! REAL(iwp) :: tempget(neq_pp1), tempput(neq_pp1)
!
! (B4) - Generate ul_pp vector
! 
  ul_pp = 0.0
  DO j = 1,nels_pp
    ul_pp(ggl_pp(:,j)) = ul_pp(ggl_pp(:,j))           &
      + utemp_pp(:)
  END DO
!
! (B5) - Generate u_pp vector - using ul_pp vectors from all
!        appropriate PEs
! 
!
! Local data:
!
  nstart = lengetsum(numpe-1)
  DO j = 1, lenget(numpe)
    u_pp(toget(j,numpe)) = u_pp(toget(j,numpe)) +        &
      ul_pp(nstart + j)
  END DO
! 
! Remote data:
! 
  DO ii = 1, numpesget
    i = pesget(ii)
    nstart = lengetsum(i-1)
    CALL MPI_ISEND (ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,             &
      MPI_COMM_WORLD,vrequest(ii),ier)
!   CALL MPI_SEND (ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,             &
!     MPI_COMM_WORLD,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) isend',ier)
  END DO
! *****
! *****    Now receive data
! *****
  DO i = 1, numpesput
    CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,     &
      vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) probe',ier)
    CALL MPI_GET_COUNT (vstatus(1,i),MPI_REAL8,recbufsize,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) get_count',ier)
    pe_number = vstatus(MPI_TAG,i)

    CALL MPI_RECV (tempput(1),lenput(pe_number),MPI_REAL8,                 &
!     MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,i),ier)
      vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) receive',ier)
    DO j = 1, lenput(pe_number)
      u_pp(toput(j,pe_number)) = u_pp(toput(j,pe_number)) + tempput(j)
    END DO
  END DO
!
! Make sure all ISENDs have completed and free up internal book-keeping of requests
!
  CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
  IF (ier /= MPI_SUCCESS) CALL mperror ('Error in scatter1 MPI_WAITALL', ier)
END SUBROUTINE scatter1
SUBROUTINE scatter2(u_pp, utemp_pp)
!
!     Performs the scatter operation: u(g) = u(g) + utemp
!
!     Final summation performed in rank order
!     Unavoidable summation still occurs on each processor
!
  IMPLICIT NONE
!
! Include file(s):
!
  INCLUDE "mpif.h"
!
!     Arguments:
!
  REAL(iwp), INTENT(INOUT) :: u_pp(neq_pp)
  REAL(iwp) :: utemp_pp(ntot,nels_pp)            ! not needed as argument(use automatic array)
  INTEGER vrequest(numpesget)
  INTEGER vstatus(MPI_STATUS_SIZE, numpesgetput)
!
!     Local Variables:
!
  INTEGER i, j, k, ii, ier, nstart, pe_number
  INTEGER ibar
  INTEGER recbufsize
  INTEGER status(MPI_STATUS_SIZE)
  INTEGER temp_status(MPI_STATUS_SIZE)
  LOGICAL lflag
  REAL(iwp), DIMENSION(0:len_pl_pp) :: ul_pp
!
! Barrier to synchronise before scattering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 23
       CALL my_barrier(numpe,ibar,details,'First barrier in gather')
!
! (B4) - Generate ul_pp vector
!  
  ul_pp = 0.0_iwp; u_pp = 0.
  DO j = 1,nels_pp
   ul_pp(ggl_pp(:,j)) = ul_pp(ggl_pp(:,j)) + utemp_pp(:,j)
  END DO
!
! (B5) - Generate u_pp vector - using ul_pp vectors from all
!        appropriate PEs
! 
! 
! Remote data:
! 
  tempput = 0.
  DO ii = 1, numpesget
    i = pesget(ii) 
  nstart = lengetsum(i-1) 
    CALL MPI_ISEND (ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,             &
      MPI_COMM_WORLD,vrequest(ii),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) isend',ier)
  END DO
! *****
! *****    Now receive data
! *****
  DO i = 1, numpesput
    CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,     &
      vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) probe',ier)
    CALL MPI_GET_COUNT (vstatus(1,i),MPI_REAL8,recbufsize,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) get_count',ier)
    pe_number = vstatus(MPI_TAG,i)
!#ifdef debug
    IF (recbufsize .NE. lenput(pe_number)) THEN
      WRITE(*,*) 'PE: ', numpe, ' Scatter error: recbufsize: ', recbufsize, ' /= lenput(pe_number): ', lenput(pe_number)
      WRITE(*,*) ' pe_number: ', pe_number, ' i: ', i
      WRITE(*,*) 'source: ', vstatus(MPI_SOURCE,i), ' tag: ', vstatus(MPI_TAG,i)
      CALL MPI_GET_COUNT (temp_status,MPI_REAL8,recbufsize,ier)
      IF (ier .NE. 0) CALL mperror('Error in (B5) get_count',ier)
      WRITE(*,*) 'Temp probe. recbufsize: ', recbufsize, ' source: ', temp_status(MPI_SOURCE),     &
        ' tag: ', temp_status(MPI_TAG)
      STOP
    ENDIF
!#endif  
    CALL MPI_RECV (tempput1(:,pe_number),lenput(pe_number),MPI_REAL8,                 &
      vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) receive',ier)
!
!   Note - need to use tempput1 to receive and store data, as must be added to
!   u_pp in rank order in the loop below. The local data is added in
!   the appropriate rank position.
!
 END DO 
!
! Make sure all ISENDs have completed and free up internal book-keeping of requests
!
  CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
  IF (ier /= MPI_SUCCESS) CALL mperror ('Error in MPI_WAITALL', ier)
   DO i=1,npes
     IF(numpe==i) THEN
       nstart = lengetsum(i-1)
       DO j = 1, lenget(i)
         u_pp(toget(j,getpes(i))) = u_pp(toget(j,getpes(i))) +           &
                                    ul_pp(nstart + j)
       END DO
     ELSE 
       DO k = 1, lenput(i)
         u_pp(toput(k,putpes(i))) = u_pp(toput(k,putpes(i))) + tempput1(k,i)
       END DO
     END IF
   END DO
!
! Barrier to synchronise after scattering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 24
       CALL my_barrier(numpe,ibar,details,'Last barrier in scatter')
END SUBROUTINE scatter2
SUBROUTINE scatter3(u_pp, utemp_pp)
!
!     Performs the scatter operation: u(g) = u(g) + utemp
!
  IMPLICIT NONE
!
! Include file(s):
!
  INCLUDE "mpif.h"  
!
!     Arguments:
!
  REAL(iwp), INTENT(INOUT) :: u_pp(neq_pp)
  REAL(iwp) :: utemp_pp(ntot,nels_pp)            ! not needed as argument(use automatic array)
  INTEGER vrequest(numpesget)
  INTEGER vstatus(MPI_STATUS_SIZE, numpesgetput)
!
!     Local Variables:
!
  INTEGER i, j, k, ii, ier, nstart, pe_number
  INTEGER ibar
  INTEGER recbufsize
  INTEGER status(MPI_STATUS_SIZE)
  INTEGER temp_status(MPI_STATUS_SIZE)
  LOGICAL lflag
  REAL(iwp), DIMENSION(0:len_pl_pp) :: ul_pp
!
! Barrier to synchronise before scattering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
!       ibar = 23
!       CALL my_barrier(numpe,ibar,details,'First barrier in gather')
!
! (B4) - Generate ul_pp vector
! 
   ul_pp = 0.0_iwp
   DO j = 1,nels_pp
     ul_pp(ggl_pp(:,j)) = ul_pp(ggl_pp(:,j)) + utemp_pp(:,j)
   END DO
! DO j = 1,nels_pp
!   DO i = 1, ntot
!     ul_pp(ggl_pp(i,j)) = ul_pp(ggl_pp(i,j))           &
!       + utemp_pp(i,j)
!   END DO
! END DO
!
! (B5) - Generate u_pp vector - using ul_pp vectors from all
!        appropriate PEs
! 
!
! Local data:
!
  nstart = lengetsum(numpe-1)
  DO j = 1, lenget(numpe)
    u_pp(toget(j,getpes(numpe))) = u_pp(toget(j,getpes(numpe))) +        &
      ul_pp(nstart + j)
  END DO
! 
! Remote data:
! 
  DO ii = 1, numpesget
    i = pesget(ii)
    nstart = lengetsum(i-1)
    CALL MPI_ISEND (ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,             &
      MPI_COMM_WORLD,vrequest(ii),ier)
!   CALL MPI_SEND (ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,             &
!     MPI_COMM_WORLD,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) isend',ier)
  END DO
! *****
! *****    Now receive data
! *****
! ibar = 500
! CALL my_barrier(numpe,ibar,details,'Scatter: Barrier before receives')
! CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,     &
!   temp_status,ier)
! IF (ier .NE. 0) CALL mperror('Error in (temp) probe',ier)
  DO i = 1, numpesput
    CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,     &
      vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) probe',ier)
    CALL MPI_GET_COUNT (vstatus(1,i),MPI_REAL8,recbufsize,ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) get_count',ier)
    pe_number = vstatus(MPI_TAG,i)
!#ifdef debug
    IF (recbufsize .NE. lenput(pe_number)) THEN
      WRITE(*,*) 'PE: ', numpe, ' Scatter error: recbufsize: ', recbufsize, ' /= lenput(pe_number): ', lenput(pe_number)
      WRITE(*,*) ' pe_number: ', pe_number, ' i: ', i
      WRITE(*,*) 'source: ', vstatus(MPI_SOURCE,i), ' tag: ', vstatus(MPI_TAG,i)
      CALL MPI_GET_COUNT (temp_status,MPI_REAL8,recbufsize,ier)
      IF (ier .NE. 0) CALL mperror('Error in (B5) get_count',ier)
      WRITE(*,*) 'Temp probe. recbufsize: ', recbufsize, ' source: ', temp_status(MPI_SOURCE),     &
        ' tag: ', temp_status(MPI_TAG)
      STOP
    ENDIF
!#endif
    CALL MPI_RECV (tempput(1),lenput(pe_number),MPI_REAL8,                 &
!     MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,i),ier)
      vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
    IF (ier .NE. 0) CALL mperror('Error in (B5) receive',ier)
!
!   Note - need to use tempput to receive data, as must be added to
!   u_pp in loop below. No problem using the same temporary buffer
!   for the receives as these are blocking receives.
!
    DO j = 1, lenput(pe_number)
      u_pp(toput(j,putpes(pe_number))) = u_pp(toput(j,putpes(pe_number))) + tempput(j)
    END DO
  END DO
!
! Make sure all ISENDs have completed and free up internal book-keeping of requests
!
  CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
  IF (ier /= MPI_SUCCESS) CALL mperror ('Error in MPI_WAITALL', ier)
!
! Barrier to synchronise after scattering data (not always required, but safer
! and removes necessity for some barriers in main program).
!
       ibar = 24
       CALL my_barrier(numpe,ibar,details,'Last barrier in scatter')
END SUBROUTINE scatter3
SUBROUTINE make_ggl(gg_pp)
!
!     Version 1a, 25-11-96, M.A. Pettipher
!     Version 1b, 11-12-96, Replaces arguments with global variables
!     Version 1c, 27-01-97, Allows more than one PE to have neq_pp2 and nels_pp2
!     Version 1d, 26-03-97, Uses vrequest and vstatus
!     Version 1e, 16-06-97, Reduces size of toget and toput by declaring second
!                           dimension as npes/4, and assigning from 1 to
!                           numpesget (or numpesput) instead of to npes. Requires
!                           trap to ensure npes/4 is big enough. Subsequently
!                           reallocates with numpesget and numpesput.
!                           Creates getpes and putpes to provide reverse
!                           indexing of pesget and pesput.
!
      IMPLICIT NONE
!
!     Include file(s):
!
      INCLUDE "mpif.h"
!
!     Generates ggl_pp and associated arrays
!
!     Arguments:
!
!     Arguments:
!
      INTEGER, INTENT(IN) :: gg_pp(ntot,nels_pp)
!
!     Local Variables:
      INTEGER i, j, k, ii, iii, ier, pe_number, bufid, count
      INTEGER position, recbufsize
      INTEGER vstatus(MPI_STATUS_SIZE, npes)        ! numpesget not yet known, so use npes
      INTEGER vrequest(npes)                        ! numpesget not yet known, so use npes
      LOGICAL lflag
      LOGICAL local                                 ! Tests whether any data is local
      INTEGER preq_pp(neq_pp1,npes)
      INTEGER rem_acc, loc_acc, sum_rem_acc, sum_numpesget
!     INTEGER :: num_neq_pp1, threshold
      REAL(iwp) :: sumtemp1, sumtemp2
!!!!! REAL(iwp) :: sumget(neq_pp1,npes_pp)
      REAL(iwp) :: rem_loc, sum_rem_loc
      LOGICAL flag, newpe
!     ***** START OF EXECUTABLE STATEMENTS *****
!
! Not sure whether this barrier is necessary
!
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      IF (ier .NE. 0) CALL mperror('Error in first barrier in make_ggl',ier)
      CALL calc_npes_pp                 ! For temporary allocation of toget
!
!     Call routine to allocate arrays for gather_scatter related operations:
!
!          WRITE(*,*) 'After flush, before allocate_gather_scatter'
          CALL allocate_gather_scatter
!          WRITE(*,*) 'After  allocate_gather_scatter'
!
          preq_pp = 0
          threshold = num_neq_pp1*neq_pp1
          DO j = 1, nels_pp
            DO i = 1,ntot
!
!             Convert index location, g, into location on specific PE.
!             preq_pp(i,j) refers to location j on ith PE.
!
!             Do not bother if index = 0 (also prevents problem
!             with assigning preq_pp(0, ))
!
              IF(gg_pp(i,j).NE.0) THEN
                IF (gg_pp(i,j) <= threshold .OR. threshold == 0) THEN
                  preq_pp(                                            &
                    MOD((gg_pp(i,j) - 1),neq_pp1) + 1,                 &
                    (gg_pp(i,j) - 1)/neq_pp1 + 1 ) = 1
                ELSE
                  preq_pp(                                            &
                    MOD((gg_pp(i,j) - threshold - 1),(neq_pp1 - 1) ) + 1,                 &
                    (gg_pp(i,j) - threshold - 1)/(neq_pp1 - 1)  + 1 + num_neq_pp1 ) = 1
                ENDIF
              ENDIF
            END DO
          END DO
! *****
! *****   (A2) - Find number of elements required from each processor (lenget(i))
! *****          Also pesget and getpes - relating actual PE number, i,  to index, ii,
! *****          of toget(:,ii).
! *****          And numpesget - number of remote PEs from which data is required.
! *****
          toget_temp = 0
          lenget = 0
          lengetsum = 0
          numpesget = 0
          getpes = 0
          ii = 0
          local = .FALSE.
          DO_find_data_outer: DO i = 1, npes
            newpe = .TRUE.
            k = 1
            DO_find_data_inner1: DO j = 1, neq_pp1
              IF(preq_pp(j,i).EQ.1) THEN
                IF (i == numpe) THEN               ! Save local until after all remote
                  local = .TRUE.
                  CYCLE DO_find_data_inner1
                END IF
                IF (newpe) ii = ii + 1
                IF (ii > npes_pp) THEN
                  WRITE(details,'(A,I6,A,I6)')                        &
                    'PE: ', numpe, ' Number of PEs to get exceeds dimension of toget: ', npes_pp
                  WRITE(details,'(A)') 'Need to increase npes_pp in main program'
                  STOP
                END IF
                newpe = .FALSE.
                toget_temp(k,ii) = j
                sumget(j,ii) = k
                k = k + 1
                lenget(i) = lenget(i) + 1
              END IF
            END DO DO_find_data_inner1
            lengetsum(i) = lengetsum(i-1) + lenget(i)
            IF (.NOT. newpe) THEN
              pesget(ii) = i
              getpes(i) = ii
              WRITE(details,'("No. of elements of PE ",I4," required by PE ",I4,": ",I8)') &
                      i, numpe, lenget(i)
            END IF
          END DO DO_find_data_outer
          numpesget = ii
          IF (local) THEN     ! For local data
            newpe = .TRUE.    ! Could avoid using newpe, but be consistent with above.
            i = numpe
            k = 1
            DO_find_data_inner2: DO j = 1, neq_pp1
              IF(preq_pp(j,i).EQ.1) THEN
                IF (newpe) ii = ii + 1
                IF (ii > npes_pp) THEN
                  WRITE(details,'(A,I6,A,I6)')                        &
                    'PE: ', numpe, ' Number of PEs to get exceeds dimension of toget: ', npes_pp
                  WRITE(details,'(A)') 'Need to increase npes_pp in main program'
                  STOP
                END IF
                newpe = .FALSE.
                toget_temp(k,ii) = j
                sumget(j,ii) = k
                k = k + 1
                lenget(i) = lenget(i) + 1
              END IF
            END DO DO_find_data_inner2
            lengetsum(i) = lengetsum(i-1) + lenget(i)
            DO iii = numpe + 1, npes
              lengetsum(iii) = lengetsum(iii) + lenget(numpe)    ! Update lengetsum for i > numpe
            END DO
            IF (.NOT. newpe) THEN
              pesget(ii) = i
              getpes(i) = ii
              WRITE(details,'("No. of elements of PE ",I4," required by PE ",I4,": ",I8)') &
                      i, numpe, lenget(i)
            END IF
          END IF

          len_pl_pp = lengetsum(npes)
          WRITE(details,*)                                              &
            'Total number of unique elements required'
          WRITE(details,*)                                              &
            'i.e. length of pl_pp required: ', len_pl_pp
          WRITE(details,'("PE: ",I4," Number of remote PEs required: ",I6)') &
           numpe, numpesget
          rem_acc = lengetsum(npes) - lenget(numpe)
          loc_acc = lenget(numpe)
          IF (loc_acc > 0) THEN
            rem_loc = rem_acc/REAL(loc_acc,iwp)
          ELSE
            rem_loc = 0
          ENDIF
          WRITE(*,'("PE: ",I4," Accesses - remote, local, remote/local: ",2I6,F8.2," From ",I6, " PEs")')  &
            numpe, rem_acc, loc_acc, rem_loc,numpesget
          CALL MPI_REDUCE(rem_loc,sum_rem_loc,1,MPI_REAL8,MPI_SUM,      &
            npes-1,MPI_COMM_WORLD,ier)
          IF(ier .NE. 0) CALL mperror('Error in REDUCE - rem_loc',ier)
          CALL MPI_REDUCE(rem_acc,sum_rem_acc,1,MPI_INTEGER,MPI_SUM,    &
            npes-1,MPI_COMM_WORLD,ier)
          IF(ier .NE. 0) CALL mperror('Error in REDUCE - rem_acc',ier)
          CALL MPI_REDUCE(numpesget,sum_numpesget,1,MPI_INTEGER,MPI_SUM,&
            npes-1, MPI_COMM_WORLD,ier)
          IF(ier .NE. 0) CALL mperror('Error in REDUCE - numpesget',ier)
          IF (numpe == npes) THEN
            WRITE(*,'("Average accesses ratio - remote/local: ", F8.2)') sum_rem_loc/REAL(npes,iwp)
            WRITE(*,'("Total remote accesses                : ", I8)') sum_rem_acc
            IF (sum_numpesget > 0)                                      &
            WRITE(*,'("Average remote accesses per PE       : ", F8.2)') sum_rem_acc/REAL(sum_numpesget,iwp)
          ENDIF

! *****
! *****   (A3) - Calculate ggl_pp
! *****
! *****   First find locations within preq_pp by summing (0 if element
! *****   not required, 1 if element required.
          DO j = 1, nels_pp
            DO i = 1,ntot
              IF (gg_pp(i,j) .EQ. 0) THEN
                ggl_pp(i,j) = 0
              ELSE
                IF (gg_pp(i,j) <= threshold .OR. threshold == 0) THEN
                  pe_number = (gg_pp(i,j) - 1)/neq_pp1 + 1
                  k = gg_pp(i,j) - (pe_number - 1)*neq_pp1
!                 k = MOD((gg_pp(i,j) - 1),neq_pp1) + 1
                ELSE
                  pe_number = num_neq_pp1 + (gg_pp(i,j) - threshold - 1)/(neq_pp1 - 1) + 1
                  k = gg_pp(i,j) - threshold - (pe_number - num_neq_pp1 - 1)*(neq_pp1 - 1)
                ENDIF
                sumtemp1 = lengetsum(pe_number - 1)
                sumtemp2 = sumget(k,getpes(pe_number))
                ggl_pp(i,j) = sumtemp1 + sumtemp2
              END IF
            END DO
          END DO
! *****
! ***** End of code to generate local index information.
! *****
! *****
! *****  (A4) - Send PE requirement request to each PE
! *****
          DO ii = 1, numpesget
            i = pesget(ii)

              CALL MPI_ISEND (toget_temp(1,ii),lenget(i),MPI_INTEGER,i-1,numpe,         &
                MPI_COMM_WORLD,vrequest(ii),ier)
              IF (ier .NE. MPI_SUCCESS) CALL mperror('Error in (A4) isend',ier)
              CALL MPI_TEST(vrequest(ii),lflag,vstatus(1,ii),ier)
!             IF(.NOT. lflag) THEN
!               WRITE (details,*) 'Send not completed:'
!               WRITE (details,*) 'source, tag:', vstatus(MPI_SOURCE,ii), vstatus(MPI_TAG,ii)
!               CALL MPI_WAIT(vrequest(ii),vstatus(1,ii),ier)
!             END IF
              WRITE(details,'("PE: ",I4," request sent to PE: ",I5)') numpe, i
!             CALL flush(details)
          END DO
! *****
! *****   (A5) - Receive PE request
! *****
! *****   Now receive corresponding data from other PEs
! *****
          CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
          IF (ier .NE. 0) CALL mperror('Error in (A5) barrier before receives',ier)
          toput_temp = 0
! *****
! *****   Test for existence of message and
! *****   receive details of elements to be sent to each PE
! *****
          count = 0
          numpesput = 0
          lenput = 0
          putpes = 0
          ii = 0
          DO                       ! Do until count limit reached
!
!           Could probably replace vstatus with status in IPROBE, GET_COUNT and RECV
!           but check carefully before changing
!
            CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,  &
              lflag,vstatus(1,numpesput+1),ier)
            IF (ier .NE. 0) CALL mperror('Error in (A5) iprobe',ier)
            IF (lflag) THEN
              count = 0
              ii = ii + 1
              numpesput = numpesput + 1         ! ii == numpesput, but prefer to use ii below
              IF (numpesput > npes_pp) THEN
                WRITE(details,'(A,I6,A,I6)')                        &
                  'PE: ', numpe, ' Number of PEs to put to exceeds dimension of toput: ', npes_pp
                STOP
              END IF           
              CALL MPI_GET_COUNT (vstatus(1,numpesput),MPI_INTEGER,recbufsize,ier)
              IF (ier .NE. 0)                                           &
                CALL mperror('Error in (A5) get_count',ier)
              pe_number = vstatus(MPI_tag,numpesput)
              lenput(pe_number) = recbufsize 
              pesput(ii) = pe_number
              putpes(pe_number) = ii
              CALL MPI_RECV (toput_temp(1,ii),lenput(pe_number),MPI_INTEGER,            &
                MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,numpesput),ier)
              IF (ier .NE. 0) CALL mperror('Error in (A5) receive',ier)
            ELSE
              count = count + 1
              IF (count .EQ. 10) EXIT          ! Allow 10 'misses'
            END IF
          END DO
!
! Make sure all ISENDs have completed and free up internal book-keeping of requests
!
          CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
          IF (ier /= MPI_SUCCESS) CALL mperror ('Error in MPI_WAITALL', ier)
!
!         Set numpesgetput to max of numpesget and numpesput
!
          numpesgetput = MAX(numpesget, numpesput)
!
!         Print more information about required communication.
!
          DO i = 1, numpesput
            WRITE(details,'("Number of elements of PE ",I4," required by PE ",I4,": ",I8)') &
              numpe, pesput(i), lenput(pesput(i))
          END DO
          WRITE(details,'("PE: ", I4," Number of PEs to send data to: ", I4)') &
            numpe, numpesput
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      IF (ier .NE. 0) CALL mperror('Error in (A5) barrier',ier)
!
! Now numpesget, numpesput, toget and toput are known, reallocate toget and toput
! Note that numpesget is number of remote PEs required, but if as would normally
! be expected, some of the data is local (on numpe), toget requires an extra element 
! for the local data. Numpesput is number to be sent to (excluding numpe).
!
      IF (local) THEN
        ALLOCATE ( toget(neq_pp1, numpesget + 1), toput(neq_pp1,numpesput) )      
        toget = toget_temp(:, 1:numpesget + 1)
      ELSE
        ALLOCATE ( toget(neq_pp1, numpesget), toput(neq_pp1,numpesput) )      
        toget = toget_temp(:, 1:numpesget)
      END IF
      toput = toput_temp(:, 1:numpesput)
      DEALLOCATE ( toget_temp, toput_temp )
!IF(numpe==1) THEN
!  DO i=1,numpesget
!      DO j=1,neq_pp
!        WRITE(12,'(a,4i5)') "toget= ",i,j,toget(j,i)
!      END DO
!  END DO
!  DO i=1,numpesput
!    DO j=1,neq_pp
!      WRITE(12,'(a,4i5)') "toput= ",i,j,toput(j,i)
!    END DO
!  END DO
!END IF
!IF(numpe==2) THEN
!  DO i=1,numpesget
!    DO j=1,neq_pp
!      WRITE(13,'(a,4i5)') "toget= ",i,j,toget(j,i)
!    END DO
!  END DO
!  DO i=1,numpesput
!    DO j=1,neq_pp
!      WRITE(13,'(a,4i5)') "toput= ",i,j,toput(j,i)
!    END DO
!  END DO
!END IF
!IF(numpe==3) THEN
!  DO i=1,numpesget
!      DO j=1,neq_pp
!        WRITE(14,'(a,4i5)') "toget= ",i,j,toget(j,i)
!      END DO
!  END DO
!  DO i=1,numpesput
!    DO j=1,neq_pp
!      WRITE(14,'(a,4i5)') "toput= ",i,j,toput(j,i)
!    END DO
!  END DO
!END IF
!IF(numpe==4) THEN
!  DO i=1,numpesget
!    DO j=1,neq_pp
!      WRITE(15,'(a,4i5)') "toget= ",i,j,toget(j,i)
!    END DO
!  END DO
!  DO i=1,numpesput
!    DO j=1,neq_pp
!      WRITE(15,'(a,4i5)') "toput= ",i,j,toput(j,i)
!    END DO
!  END DO
!END IF
!
! The following barrier is probably necessary, but not checked thoroughly
!
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      IF (ier .NE. 0) CALL mperror('Error in last barrier in make_ggl',ier)
END SUBROUTINE make_ggl
END MODULE gather_scatter6
