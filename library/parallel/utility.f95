MODULE utility

!
! Contains some general routines for message passing and timing operations
!

  USE global_variables1 
  USE precision
  USE mp_module
  IMPLICIT NONE


CONTAINS
      SUBROUTINE shutdown()
      IMPLICIT NONE
      CALL MPI_FINALIZE(ier)
      RETURN
      END SUBROUTINE shutdown


      SUBROUTINE my_barrier(numpe,ibar,channel,chstr)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: numpe, ibar, channel
      CHARACTER (LEN=*), INTENT(IN) :: chstr
      INTEGER :: ier



      INCLUDE "mpif.h"





      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      IF (ier .NE. 0) CALL mperror(chstr,ier)


      IF (ier .NE. 0) CALL mperror(chstr,ier)
      END SUBROUTINE my_barrier

!--------------------------------------------------------------------------

      SUBROUTINE find_pe_procs(numpe,numprocs)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: numpe, numprocs




      INTEGER :: ier
      INCLUDE "mpif.h"
      CALL MPI_INIT(ier)
!      CALL mpi_nail()
      IF (ier .NE. 0) CALL mperror('Error in MPI_INIT',ier)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,numpe,ier)
      IF (ier .NE. 0) CALL mperror('Error in MPI_COMM_RANK',ier)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ier)
      IF (ier .NE. 0) CALL mperror('Error in MPI_COMM_SIZE',ier)

      numpe = numpe + 1                   ! From 1 to npes
      END SUBROUTINE find_pe_procs

!--------------------------------------------------------------------------

      SUBROUTINE mperror(cherror, errcode)
      IMPLICIT NONE
      CHARACTER*(*) cherror
      INTEGER errcode


!     INTEGER MPI_Max_error_string
!     CHARACTER*(MPI_Max_error_string) errstring


      WRITE(*,'(A,A,I5)') cherror, ' errcode = ', errcode






      STOP
      END SUBROUTINE mperror

!--------------------------------------------------------------------------

       FUNCTION elap_time()
       REAL(iwp) elap_time
 
 
       INTEGER count, count_rate
       CALL SYSTEM_CLOCK(count,count_rate)
       elap_time = REAL(count,iwp)/REAL(count_rate,iwp)
 
 
       END FUNCTION elap_time

!--------------------------------------------------------------------------
 
 !     FUNCTION cpu_time()
 !     REAL(iwp) cpu_time, cpu
 
 
 !      cpu_time = 0.0
 
 
 !      END FUNCTION cpu_time
 
!--------------------------------------------------------------------------

        FUNCTION reduce(neq_temp) RESULT(neq)
        ! maximum neq on any processor
        IMPLICIT NONE; INTEGER, INTENT(IN)::neq_temp; INTEGER::neq
        CALL MPI_ALLREDUCE(neq_temp,neq,                                 &
                           1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
        END FUNCTION reduce

        FUNCTION dot_product_p(vectora,vectorb)

        !calculates dot_product of distributed vectors
        !different summation order to serial fortran intrinsic dot_product
 
        REAL(iwp), INTENT(IN)::vectora(:),vectorb(:)
        REAL(iwp)::local_product,dot_product_p
        
        local_product=dot_product(vectora,vectorb)
        call MPI_ALLREDUCE(local_product,dot_product_p,1,MPI_REAL8,      &
                           MPI_SUM,MPI_COMM_WORLD,ier)

        END FUNCTION dot_product_p

!-------------------------------------------------------------------------
  
        FUNCTION norm_p(vectora)

        !calculates 12norm of distributed vector
        !different summation order to serial version

        REAL(iwp), INTENT(IN)::vectora(:)
        REAL(iwp)::local_product,global_product,norm_p
 
        local_product=dot_product(vectora,vectora)
        call MPI_ALLREDUCE(local_product,global_product,1,MPI_REAL8,     &
                           MPI_SUM,MPI_COMM_WORLD,ier)
        norm_p=sqrt(global_product)
  
        END FUNCTION norm_p

!------------------------------------------------------------------------

        FUNCTION sum_p(vectora)

        !sums a vector distributed over processors

        REAL(iwp), INTENT(IN)::vectora(:)
        REAL(iwp)::local_sum, global_sum , sum_p

        local_sum = SUM(vectora)
        CALL MPI_ALLREDUCE(local_sum,global_sum,1,MPI_REAL8,MPI_SUM,    &
                           MPI_COMM_WORLD,ier)
        sum_p = global_sum

        END FUNCTION sum_p

!---------------------------------------------------------------------------

        FUNCTION maxval_p(vectora)

        !maximum value of real distributed vector  

        REAL(iwp), INTENT(IN)::vectora(:)
        REAL(iwp)::local_max, global_max , maxval_p

        local_max = MAXVAL(vectora)
        CALL MPI_ALLREDUCE(local_max,global_max,1,MPI_REAL8,MPI_SUM,    &
                           MPI_COMM_WORLD,ier)
        maxval_p = global_max

        END FUNCTION maxval_p

!---------------------------------------------------------------------------

      SUBROUTINE readbf(nf,nn,nodof,nrb)
    
      ! 'set' input of nf data
    
      IMPLICIT NONE
      INTEGER nn,nodof,nrb,nf(nn,nodof),nfd(nodof)
      INTEGER i,j,l,n,if,it,is
      nf=1
      DO l=1,nrb
         READ(10,*)if,it,is,nfd
         DO i=if,min(nn,it),is 
           nf(i,:) = nf(i,:) * nfd(:)
         END DO
      END DO
      n=0
    
      DO i=1,nn 
        DO j=1,nodof

          l=nf(i,j)
          if(nf(i,j) .NE. 0) THEN
            n=n+1
            nf(i,j)=n
          END IF
        END DO
      END DO
      RETURN
    
      END SUBROUTINE readbf



SUBROUTINE reindex_fixed_nodes(ieq_start,no,no_local_temp,  &
                               num_no,no_index_start)

! Local equation indexing

IMPLICIT NONE

integer, INTENT(IN)     :: ieq_start, no(:)
integer, INTENT(INOUT)  :: no_local_temp(:)
integer, INTENT(OUT)    :: num_no, no_index_start
integer                 :: fixed_nodes, i

fixed_nodes = ubound(no,1)

  no_index_start = 0
  num_no = 0
   do i=1,fixed_nodes
     if (no(i) < ieq_start) then
       cycle
     else if (no(i) >= ieq_start + neq_pp) then
       exit
     else if (no_index_start==0) then
       no_index_start = i
       no_local_temp(1) = no(i)
       num_no = 1
     else
       num_no = num_no + 1
       no_local_temp(num_no) = no(i)
     end if
   end do
 
END SUBROUTINE reindex_fixed_nodes

SUBROUTINE glob_sum_sq(neq,ieq_start,vectora_pp,vectorb_pp,answer)

!------------------------------------------------------------------------------
!---- Performs a sum of squares for real vectors that are on different
!---- processors ie performs sum of (vectora(i)*vectorb(i))
!------------------------------------------------------------------------------

implicit none

integer               :: neq,ieq_start
real(iwp),intent(out) :: answer
real(iwp),intent(in)  :: vectora_pp(:),vectorb_pp(:)
real(iwp),allocatable :: vectora_tmp(:), vectorb_tmp(:), vectora(:), vectorb(:)


allocate(vectora_tmp(neq),vectorb_tmp(neq),vectora(neq),vectorb(neq))

vectora_tmp = 0. ; vectorb_tmp = 0. ; vectora = .0 ; vectorb = 0.

vectora_tmp(ieq_start:ieq_start+neq_pp-1) = vectora_pp
vectorb_tmp(ieq_start:ieq_start+neq_pp-1) = vectorb_pp

call MPI_ALLREDUCE(vectora_tmp,vectora,neq,MPI_REAL8,MPI_SUM,    &
                   MPI_COMM_WORLD,ier)
call MPI_Barrier(MPI_COMM_WORLD,ier)
call MPI_ALLREDUCE(vectorb_tmp,vectorb,neq,MPI_REAL8,MPI_SUM,    &
                   MPI_COMM_WORLD,ier)                   
call MPI_Barrier(MPI_COMM_WORLD,ier)

deallocate(vectora_tmp, vectorb_tmp)

answer=dot_product(vectora,vectorb)

deallocate(vectora,vectorb)

return

END SUBROUTINE glob_sum_sq

SUBROUTINE construct_full(x,ieq_start,xfull,xoldfull)

!
! constructs a full (serial) vector from a distributed (parallelised) vector
!

implicit none

integer::neq, neq_pp
integer, intent(in)::ieq_start
real(iwp),intent(in):: x(:)
real(iwp),intent(inout):: xfull(:),xoldfull(:)
real(iwp),allocatable:: xfulltmp(:)

neq=ubound(xfull,1); neq_pp=ubound(x,1)

allocate(xfulltmp(SIZE(xfull)))
xfulltmp=0.
xoldfull = xfull; xfull = 0.
xfulltmp(ieq_start+1:ieq_start+neq_pp) = x
call MPI_ALLREDUCE(xfulltmp(2:neq+1),xfull(2:neq+1),neq,MPI_REAL8,MPI_SUM,    &
                   MPI_COMM_WORLD,ier)
call MPI_Barrier(MPI_COMM_WORLD,ier)
deallocate(xfulltmp)

return

END SUBROUTINE construct_full
 
SUBROUTINE constructv(xfull,ieq_start,x,numpe)

!
! constructs a full (serial) vector from a distributed (parallelised) vector
! then returns the result to the processor that called
!  

implicit none

integer::neq, neq_pp
integer, intent(in)::ieq_start,numpe
real(iwp),intent(in):: x(:)
real(iwp),intent(inout):: xfull(:)
real(iwp),allocatable:: xfulltmp(:)

neq=ubound(xfull,1); neq_pp=ubound(x,1)
allocate(xfulltmp(SIZE(xfull)))
xfulltmp=0.
xfull = 0.
xfulltmp(ieq_start+1:ieq_start+neq_pp) = x
call MPI_ALLREDUCE(xfulltmp(2:neq+1),xfull(2:neq+1),neq,MPI_REAL8,MPI_SUM,    &
                   MPI_COMM_WORLD,ier)
call MPI_Barrier(MPI_COMM_WORLD,ier)
deallocate(xfulltmp)

return

END SUBROUTINE constructv


SUBROUTINE checon_par(loads,oldlds,tol,converged,neq_pp)

!
! parallel version of checon
! sets converged to .false. if relative change in loads and
! oldlds is greater than tol and updates oldlds
!

IMPLICIT NONE

real(iwp)::maxloads_pp,maxdiff_pp,maxloads,maxdiff
real(iwp),intent(in)::loads(:),tol;real(iwp),intent(inout)::oldlds(:)
logical,intent(out)::converged
integer::i,neq_pp 
  maxloads_pp=0.; maxdiff_pp=0.  ; converged = .true.

  do i=1,neq_pp
    if(abs(loads(i))>maxloads_pp) maxloads_pp=abs(loads(i)) 
  end do
  CALL MPI_ALLREDUCE(maxloads_pp,maxloads,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ier)
  
  do i=1,neq_pp
    if(abs((loads(i)-oldlds(i))/maxloads)>maxdiff_pp) then
      maxdiff_pp=abs((loads(i)-oldlds(i))/maxloads)
    end if 
  end do

  CALL MPI_ALLREDUCE(maxdiff_pp,maxdiff,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ier)
   
  if (maxdiff>tol) converged=.false. 
  oldlds=loads
 return

!
! not sure if we need barriers or something, need to check if mpi_allreduce
! puts in a barrier automatically
!
! initially keep the large vector loads/oldloads for other steps to check

END SUBROUTINE checon_par

SUBROUTINE bcast_inputdata_p121 (numpe,npes,nels,nxe,nze,nip,            &
               aa,bb,cc,e,v,tol,limit)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,limit  
real(iwp), INTENT(INOUT)   ::aa,bb,cc,e,v,tol
!
! Assign temporary buffer for broadcast data
!

bufsizer=5*ilength + 6*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                    &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,                     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
end if

END SUBROUTINE bcast_inputdata_p121     

SUBROUTINE bcast_inputdata_p122(numpe,npes,phi,c,psi,e,v,cons, nels,nxe,nze,&
                                nip,aa,bb,cc,incs,plasits,cjits,plastol,cjtol)   
IMPLICIT NONE 
Integer , INTENT(INOUT):: numpe,npes,nels,nxe,nze,nip,incs,plasits,cjits
Real(iwp), INTENT(INOUT) ::    phi,c,psi,e,v,cons,aa,bb,cc,plastol,cjtol

bufsizer = 7 *ilength + 11 *rlength 
CALL MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)
bufdecl=bufsizer/4   ; ALLOCATE(tempbuf(bufdecl))   
IF(numpe==npes) THEN
position = 0
!--------- INTEGERS -------------------

CALL MPI_PACK(nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK(nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK(nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK(nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK(incs,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (plasits,1,MPI_INTEGER,tempbuf,bufsizer,position,             &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (cjits,1,MPI_INTEGER,tempbuf,bufsizer,position,               &
              MPI_COMM_WORLD,ier) 

!-----REALS----------------------------

CALL MPI_PACK (phi,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (c,1,MPI_REAL8,tempbuf,bufsizer,position,                     &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (psi,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                     &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                     &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (cons,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                    &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                    &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                    &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (plastol,1,MPI_REAL8,tempbuf,bufsizer,position,               &
              MPI_COMM_WORLD,ier)
CALL MPI_PACK (cjtol,1,MPI_REAL8,tempbuf,bufsizer,position,                 &
              MPI_COMM_WORLD,ier) 
END IF

CALL MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

IF(numpe/=npes) THEN  

position = 0

!--------- INTEGERS --------------------

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,            &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,              &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,              &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,incs,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,plasits,1,MPI_INTEGER,            &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjits,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !

!---- REALS -----------------------

CALL MPI_UNPACK (tempbuf,bufsizer,position,phi,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,c,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,psi,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cons,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,plastol,1,MPI_REAL8,            &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjtol,1,MPI_REAL8,              &
     &      MPI_COMM_WORLD,ier) !
END IF

END SUBROUTINE bcast_inputdata_p122 
 

SUBROUTINE bcast_inputdata_p123 (numpe,npes,nels,nxe,nze,nip,              &
              aa,bb,cc,kx,ky,kz,tol,limit,loaded_freedoms,fixed_freedoms)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,limit,               &
                        loaded_freedoms,fixed_freedoms
real(iwp), INTENT(INOUT)   ::aa,bb,cc,tol,kx,ky,kz
!
! Assign temporary buffer for broadcast data
!

bufsizer=7*ilength + 7*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (loaded_freedoms,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (fixed_freedoms,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kx,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (ky,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kz,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,loaded_freedoms,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,fixed_freedoms,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kx,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,ky,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kz,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p123


SUBROUTINE bcast_inputdata_p124 (numpe,npes,nels,nxe,nze,nip,             &
              aa,bb,cc,kx,ky,kz,dtim,nstep,theta,npri,             &
              tol,limit,val0)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,nstep,npri,limit               
real(iwp), INTENT(INOUT)   ::aa,bb,cc,tol,kx,ky,kz,dtim,theta,val0
!
! Assign temporary buffer for broadcast data
!

bufsizer=7*ilength + 10*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (npri,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kx,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (ky,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kz,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (dtim,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (theta,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (val0,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,npri,1,MPI_INTEGER,               & 
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kx,1,MPI_REAL8,              &  
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,ky,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kz,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,dtim,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,theta,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,val0,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p124

SUBROUTINE bcast_inputdata_p125 (numpe,npes,nels,nxe,nze,nip,            &
              aa,bb,cc,kx,ky,kz,dtim,nstep,npri,val0)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,nstep,npri         
real(iwp), INTENT(INOUT)   ::aa,bb,cc,kx,ky,kz,dtim,val0
!
! Assign temporary buffer for broadcast data
!

bufsizer=6*ilength + 8*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (npri,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kx,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (ky,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kz,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (dtim,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (val0,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,npri,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kx,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,ky,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kz,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,dtim,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,val0,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p125

SUBROUTINE bcast_inputdata_p95 (numpe,npes,nels,nxe,nze, aa,bb,cc,       &
                    nr,nip,permx,permy,permz,e,v,dtim,nstep,theta,       &
                       cjits,cjtol,loaded_freedoms)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nr,nip,nstep,cjits,      &
                        loaded_freedoms
real(iwp), INTENT(INOUT)   ::aa,bb,cc,permx,permy,permz,e,v,dtim,theta,cjtol
!
! Assign temporary buffer for broadcast data
!

bufsizer=8*ilength + 11*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nr,1,MPI_INTEGER,tempbuf,bufsizer,position,                    &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjits,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (loaded_freedoms,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (permx,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (permy,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (permz,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (dtim,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (theta,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjtol,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nr,1,MPI_INTEGER,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjits,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,loaded_freedoms,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,permx,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,permy,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,permz,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,dtim,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,theta,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjtol,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p95

SUBROUTINE bcast_inputdata_p127 (numpe,npes,nels,nxe,nze, aa,bb,cc,       &
                        nip,kx,ky,kz,e,v,dtim,nstep,theta,cjits,cjtol)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,nstep,cjits  
real(iwp), INTENT(INOUT)   ::aa,bb,cc,kx,ky,kz,e,v,dtim,theta,cjtol
!
! Assign temporary buffer for broadcast data
!

bufsizer=6*ilength + 11*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjits,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kx,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (ky,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kz,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (dtim,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (theta,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjtol,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjits,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kx,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,ky,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kz,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,dtim,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,theta,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjtol,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p127

SUBROUTINE bcast_inputdata_p96 (numpe,npes,nels,nxe,nze, nr,nip,aa,bb,cc,   &
                                fixed_nodes,visc,rho,tol,limit,cjtol,cjits, &
                                penalty,x0,ell,kappa,                       &
                                nrb,nrb1,nrb2,nrb3,nrb4,nrb5,nrb6)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nr,nip,fixed_nodes,      &
                        limit,cjits,ell,nrb,nrb1,nrb2,nrb3,nrb4,nrb5,nrb6      
real(iwp), INTENT(INOUT)   ::aa,bb,cc,visc,rho,tol,cjtol,penalty,x0,kappa
!
! Assign temporary buffer for broadcast data
!

bufsizer=16*ilength + 10*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,           &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,            &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,            &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nr,1,MPI_INTEGER,tempbuf,bufsizer,position,             &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,            &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (fixed_nodes,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjits,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (ell,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nrb,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nrb1,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nrb2,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nrb3,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nrb4,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nrb5,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nrb6,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (visc,1,MPI_REAL8,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjtol,1,MPI_REAL8,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (penalty,1,MPI_REAL8,tempbuf,bufsizer,position,      &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (x0,1,MPI_REAL8,tempbuf,bufsizer,position,           &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kappa,1,MPI_REAL8,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nr,1,MPI_INTEGER,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,fixed_nodes,1,MPI_INTEGER, &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjits,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,ell,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nrb,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nrb1,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nrb2,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nrb3,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nrb4,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nrb5,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nrb6,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,visc,1,MPI_REAL8,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,        &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjtol,1,MPI_REAL8,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,penalty,1,MPI_REAL8,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,x0,1,MPI_REAL8,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kappa,1,MPI_REAL8,        &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p96

SUBROUTINE bcast_inputdata_p126(numpe,npes,nels,nxe,nze,nip,aa,bb,cc,   &
                       visc,rho,tol,limit,cjtol,cjits,penalty,x0,ell,kappa)                       
                                

IMPLICIT NONE


integer,INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,limit,cjits,ell      
real(iwp), INTENT(INOUT)   ::aa,bb,cc,visc,rho,tol,cjtol,penalty,x0,kappa
!
! Assign temporary buffer for broadcast data
!

bufsizer=7*ilength + 10*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,           &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,            &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,            &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,            &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjits,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (ell,1,MPI_INTEGER,tempbuf,bufsizer,position,     &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (visc,1,MPI_REAL8,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjtol,1,MPI_REAL8,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (penalty,1,MPI_REAL8,tempbuf,bufsizer,position,      &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (x0,1,MPI_REAL8,tempbuf,bufsizer,position,           &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kappa,1,MPI_REAL8,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjits,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,ell,1,MPI_INTEGER,  &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,         &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,visc,1,MPI_REAL8,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,        &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjtol,1,MPI_REAL8,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,penalty,1,MPI_REAL8,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,x0,1,MPI_REAL8,          &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kappa,1,MPI_REAL8,        &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p126

SUBROUTINE bcast_inputdata_p104 (numpe,npes,nels,nye,nr,nip,aa,bb,   &
                                 rho,e,v,nmodes,el,er)

IMPLICIT NONE


integer,INTENT(INOUT)::numpe,npes,nels,nye,nr,nip,nmodes
real(iwp), INTENT(INOUT)   ::aa,bb,rho,e,v,el,er
!
! Assign temporary buffer for broadcast data
!

bufsizer=5*ilength + 7*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nye,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nr,1,MPI_INTEGER,tempbuf,bufsizer,position,                    &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nmodes,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (el,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (er,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nye,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nr,1,MPI_INTEGER,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nmodes,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,el,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,er,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p104


SUBROUTINE bcast_inputdata_p128 (numpe,npes,nels,nxe,nze,nip,aa,bb,cc,  &
                      rho,e,v,nmodes,el,er,lalfa,leig,lx,lz,acc)

IMPLICIT NONE


integer,INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,nmodes,lalfa,leig,lx,lz 
real(iwp), INTENT(INOUT)   ::aa,bb,cc,rho,e,v,el,er,acc
!
! Assign temporary buffer for broadcast data
!

bufsizer=9*ilength + 9*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nmodes,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (lalfa,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (leig,1,MPI_INTEGER,tempbuf,bufsizer,position,                    &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (lx,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (lz,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (el,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (er,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (acc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nmodes,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,lalfa,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,leig,1,MPI_INTEGER,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,lx,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,lz,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,el,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,er,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,acc,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p128

SUBROUTINE bcast_inputdata_p117 (numpe,npes,nels,nxe,nze,nr,nip,            &
              aa,bb,cc,rho,e,v,alpha1,beta1,nstep,theta,npri,               &
                       omega,tol,limit)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nr,nip,nstep,npri,limit           
real(iwp), INTENT(INOUT)   ::aa,bb,cc,rho,e,v,alpha1,beta1,theta,omega,tol
!
! Assign temporary buffer for broadcast data
!

bufsizer=8*ilength + 11*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nr,1,MPI_INTEGER,tempbuf,bufsizer,position,                    &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (npri,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (theta,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (alpha1,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (beta1,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (omega,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nr,1,MPI_INTEGER,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,npri,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,theta,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,alpha1,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,beta1,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,omega,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p117

SUBROUTINE bcast_inputdata_p129 (numpe,npes,nels,nxe,nze,nip,            &
              aa,bb,cc,rho,e,v,alpha1,beta1,nstep,theta,npri,               &
                       omega,tol,limit)

IMPLICIT NONE


integer, INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,nstep,npri,limit           
real(iwp), INTENT(INOUT)   ::aa,bb,cc,rho,e,v,alpha1,beta1,theta,omega,tol
!
! Assign temporary buffer for broadcast data
!

bufsizer=7*ilength + 11*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (npri,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (theta,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (alpha1,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (beta1,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (omega,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,npri,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,theta,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,alpha1,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,beta1,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,omega,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p129

SUBROUTINE bcast_inputdata_p118 (numpe,npes,nels,nxe,nze,nr,nip,loaded_nodes,&
              aa,bb,cc,rho,e,v,sbary,pload,nstep,dtim,npri)

IMPLICIT NONE


integer,INTENT(INOUT)::numpe,npes,nels,nxe,nze,nr,nip,loaded_nodes,nstep,npri
real(iwp), INTENT(INOUT)   ::aa,bb,cc,rho,e,v,sbary,pload,dtim
!
! Assign temporary buffer for broadcast data
!

bufsizer=8*ilength + 9*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nr,1,MPI_INTEGER,tempbuf,bufsizer,position,                    &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (loaded_nodes,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (npri,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (sbary,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (dtim,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (pload,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nr,1,MPI_INTEGER,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,loaded_nodes,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,npri,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,sbary,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,dtim,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,pload,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p118

SUBROUTINE bcast_inputdata_p1210 (numpe,npes,nels,nxe,nze,nip,loaded_nodes,&
              aa,bb,cc,rho,e,v,sbary,pload,nstep,dtim,npri)

IMPLICIT NONE


integer,INTENT(INOUT)::numpe,npes,nels,nxe,nze,nip,loaded_nodes,nstep,npri
real(iwp), INTENT(INOUT)   ::aa,bb,cc,rho,e,v,sbary,pload,dtim
!
! Assign temporary buffer for broadcast data
!

bufsizer=7*ilength + 9*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (loaded_nodes,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (npri,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (sbary,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (rho,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (dtim,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (pload,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,loaded_nodes,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,npri,1,MPI_INTEGER,      &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,sbary,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,rho,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,dtim,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,pload,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p1210  

END MODULE utility
