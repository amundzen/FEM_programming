MODULE timing
  USE precision
  IMPLICIT NONE
!
!     Real variables for timing:
!
  REAL(iwp) :: time(20) = 0.0_iwp, timc(20) = 0.0_iwp
  REAL(iwp) :: timest(20), timcst(20), timeend(20), timcend(20)
  REAL(iwp) :: maxtime(20) = 0.0_iwp, maxtimc(20) = 0.0_iwp, &
               mintime(20) = 999999.0_iwp, mintimc(20) = 999999.0_iwp

  REAL(IWP), PARAMETER :: c2sec = 1.0_iwp
  REAL(iwp) :: flops, flops1, flops2, mflops, mflops1, mflops2
  REAL(iwp) :: start_time,end_time
!
!     Variables used for summary file
!
  INTEGER, PARAMETER :: summary = 20
  CHARACTER (LEN=8)  :: chdate
  CHARACTER (LEN=10) :: chtime
  CHARACTER (LEN=47) :: summary_file
  CHARACTER (LEN=6)  :: chnpes
  CHARACTER (LEN=6)  :: chnels

CONTAINS

SUBROUTINE getfilename(nels,npes,summary_file)

!
! Create file name for summary output file
! Include number of PEs, data size and date in file name
!

IMPLICIT NONE
 
INTEGER, INTENT(IN) ::    nels, npes
CHARACTER, INTENT(OUT) :: summary_file*47

 CALL DATE_AND_TIME(chdate,chtime)
  SELECT CASE (nels)
    CASE(216)
      chnels = 'c06'
    CASE(1000)
      chnels = 'c10'
    CASE (8000)
      chnels = 'c20'
    CASE (27000)
      chnels = 'c30'
    CASE DEFAULT
      chnels = 'nels'
!     WRITE (*,*) 'Incorrect number of elements'
!     STOP
  END SELECT
  WRITE (chnpes,'(I6)') npes
  WRITE (*,'(A)' ) chnels
  WRITE(*,'(5A)') chnpes, ' - ', ADJUSTL(chnpes),' - ', TRIM(ADJUSTL(chnpes))
  summary_file = 'par903zsmv_checon_par.summary_'//TRIM(ADJUSTL(chnpes))//'_'//TRIM(chnels)//  &
   '_'//chdate(7:8)//'-'//chdate(5:6)//'_'//chtime(1:2)//':'//chtime(3:4)//  &
   ':'//chtime(5:6)//'.'//chtime(8:10)

END SUBROUTINE getfilename

SUBROUTINE print_timing_info(numpe,npes,iters,cjiters,ell,nels,nels_pp)

IMPLICIT NONE

INCLUDE "mpif.h"

integer::ier
integer, INTENT(IN)::numpe,npes,iters,cjiters,ell,nels,nels_pp

!
! Find maximum times:
!

  CALL MPI_ALLREDUCE(time(1),maxtime(1),20,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ier)
  CALL MPI_ALLREDUCE(timc(1),maxtimc(1),20,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ier)
  
  CALL MPI_ALLREDUCE(time(1),mintime(1),20,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ier)
  CALL MPI_ALLREDUCE(timc(1),mintimc(1),20,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ier)


 !IF (numpe == npes) THEN
 !  WRITE(summary,'(//"Maximum Timing Details:"/)')
 !  WRITE(summary,'("PE: ",I4," Total time     : ", E16.4)') numpe, maxtime(1)
 !  WRITE(summary,'("PE: ",I4," Elements3      : ", E16.4)') numpe, maxtime(2)
 !  WRITE(summary,'("PE: ",I4," Elements4      : ", E16.4)') numpe, maxtime(3)
 !  WRITE(summary,'("PE: ",I4," Elements5      : ", E16.4)') numpe, maxtime(4)
 !  flops = REAL( iters*(1+2*cjiters*ell)*nels_pp*( (2*8*60)+(2*60*8)+(2*20*20)+(2*20*20)+(2*20*20) + 68), iwp)
 !  mflops = 1.0E-6*flops/( maxtime(2) + maxtime(3) + maxtime(4) )
 !  WRITE(summary,'("PE: ",I4," Elements345 Mflop/s: ", E16.4)') numpe, mflops
 !END IF

IF (numpe==npes) THEN
  WRITE(summary,'("SUMMARY PERFORMANCE OF PROGRAM PAR903ZSMVO"/              &
                & "COMPUTER   : FERMAT "/"ELEMENTS   : ",I6,/"PROCESSORS : " &
                & ,a,"       DATE: ",a,"-",a,"-",a,"  TIME: ",a,"-",a/)')    &
                & nels, chnpes, chdate(7:8), chdate(5:6), chdate(3:4),       &
                & chtime(1:2), chtime(3:4)


  WRITE(summary,'("PROGRAM TOTALS:"/"~~~~~~~~~~~~~~~"/)')
  WRITE(summary,'(" Program time       Main loop - BiCG    BiCGSTAB            &
       &Other operations"/)')
  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              & "    user: ",E10.4)') maxtime(1), maxtime(2)-maxtime(3),       &
                maxtime(3),maxtime(1) - maxtime(2)
  WRITE(summary,'("%user: ","       -  ","  %user: ",f10.2,"   %user: ",f10.2, &
              & "   %user: ",f10.2)') ((maxtime(2)-maxtime(3))/maxtime(1))*100,&
              (maxtime(3)/maxtime(1))*100,                                     &
              ((maxtime(1) - maxtime(2))/maxtime(1))*100                 
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              & "     cpu: ",E10.4)') maxtimc(1), maxtimc(2)-maxtimc(3),       &
              maxtimc(3), maxtimc(1) - maxtimc(2)  
  WRITE(summary,'(" %cpu: ","       -  ","   %cpu: ",f10.2,                    &
             &"    %cpu: ",f10.2,"    %cpu: ",f10.2/)')                        &
             ((maxtimc(2)-maxtimc(3))/maxtimc(1))*100,                         &
              (maxtimc(3)/maxtimc(1))*100,                                     &
              ((maxtimc(1) - maxtimc(2))/maxtimc(1))*100                 

  WRITE(summary,'("BREAKDOWN OF BiCGSTAB:"/"~~~~~~~~~~~~~~~~~~~~~~"/)')

  WRITE(summary,'(" Initialisation     Elements_45         Other operations    &
       &total (c.f. above)"//" MAXIMUMS")')
       
  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              & "    user: ",E10.4)')                                          &
       maxtime(5), maxtime(6)+maxtime(7)+maxtime(8)+maxtime(9),                &
       maxtime(3)-maxtime(5)-maxtime(6)-maxtime(7)-maxtime(8)-maxtime(9),      &
       maxtime(3)  
  WRITE(summary,'("%user: ",f10.2,"  %user: ",f10.2,"   %user: ",f10.2,        &
              & "   %user: ","       -  ")')                                   &
     (maxtime(5)/maxtime(3))*100,                                              &
     ((maxtime(6)+maxtime(7)+maxtime(8)+maxtime(9))/maxtime(3))*100,           &
     ((maxtime(3)-maxtime(5)-maxtime(6)-maxtime(7)-maxtime(8)                  &
       -maxtime(9))/maxtime(3))*100 
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              & "     cpu: ",E10.4)')                                          &
       maxtimc(5), maxtimc(6)+maxtimc(7)+maxtimc(8)+maxtimc(9),                &
       maxtimc(3)-maxtimc(5)-maxtimc(6)-maxtimc(7)-maxtimc(8)-maxtimc(9),      &
       maxtimc(3)  
  WRITE(summary,'(" %cpu: ",f10.2,"   %cpu: ",f10.2,                           &
               &"    %cpu: ",f10.2,"    %cpu: ","       -  "/)')               &
     (maxtimc(5)/maxtimc(3))*100,                                              &
     ((maxtimc(6)+maxtimc(7)+maxtimc(8)+maxtimc(9))/maxtimc(3))*100,           &
     ((maxtimc(3)-maxtimc(5)-maxtimc(6)-maxtimc(7)-maxtimc(8)                  &
       -maxtimc(9))/maxtimc(3))*100 
 
 WRITE(summary,'(" MINIMUMS")')                                          
  
 WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,         &
              & "    user: ",E10.4)')                                          &
       mintime(5),mintime(6)+mintime(7)+mintime(8)+mintime(9),                 &
       mintime(3)-mintime(5)-mintime(6)-mintime(7)-mintime(8)-mintime(9),      &
       mintime(3)  
  WRITE(summary,'("%user: ",f10.2,"  %user: ",f10.2,"   %user: ",f10.2,        &
              & "   %user: ","       -  ")')                                   &
     (mintime(5)/mintime(3))*100,                                              &
     ((mintime(6)+mintime(7)+mintime(8)+mintime(9))/mintime(3))*100,           &
     ((mintime(3)-mintime(5)-mintime(6)-mintime(7)-mintime(8)                  &
       -mintime(9))/mintime(3))*100 
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              & "     cpu: ",E10.4)')                                          &
       mintimc(5), mintimc(6)+mintimc(7)+mintimc(8)+mintimc(9),                &
       mintimc(3)-mintimc(5)-mintimc(6)-mintimc(7)-mintimc(8)-mintimc(9),      &
       mintimc(3)  
  WRITE(summary,'(" %cpu: ",f10.2,"   %cpu: ",f10.2,                           &
               &"    %cpu: ",f10.2,"    %cpu: ","       -  "/)')               &
      (mintimc(5)/mintimc(3))*100,                                             & 
     ((mintimc(6)+mintimc(7)+mintimc(8)+mintimc(9))/mintimc(3))*100,           &
     ((mintimc(3)-mintimc(5)-mintimc(6)-mintimc(7)-mintimc(8)                  &
       -mintimc(9))/mintimc(3))*100 
 
  WRITE(summary,'("BREAKDOWN OF ELEMENTS_4:"/"~~~~~~~~~~~~~~~~~~~~~~~~"/)')
 
  WRITE(summary,'(" Gather             Matmul              Scatter             &
       &total (c.f. below)"//" MAXIMUMS")')
 
  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              &  "    user: ",E10.4)') maxtime(6), maxtime(7), maxtime(8),     &
                  maxtime(6)+maxtime(7)+maxtime(8)
  WRITE(summary,'("%user: ",f10.2,"  %user: ",f10.2,"   %user: ",f10.2,        &
              &  "   %user: ","       -  ")')                                  &
              (maxtime(6)/(maxtime(6)+maxtime(7)+maxtime(8)))*100,             &
              (maxtime(7)/(maxtime(6)+maxtime(7)+maxtime(8)))*100,             &
              (maxtime(8)/(maxtime(6)+maxtime(7)+maxtime(8)))*100
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              &  "     cpu:  ",E10.4)') maxtimc(6), maxtimc(7), maxtimc(8),    &
                  maxtimc(6)+maxtimc(7)+maxtimc(8)
  WRITE(summary,'(" %cpu: ",f10.2,"   %cpu: ",f10.2,                           &
             &"    %cpu:  ",f10.2,"   %cpu: ","       -  "/)')                 &
              (maxtimc(6)/(maxtimc(6)+maxtimc(7)+maxtimc(8)))*100,             &
              (maxtimc(7)/(maxtimc(6)+maxtimc(7)+maxtimc(8)))*100,             &
              (maxtimc(8)/(maxtimc(6)+maxtimc(7)+maxtimc(8)))*100

  WRITE(summary,'(" MINIMUMS")')

  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              &  "    user: ",E10.4)') mintime(6), mintime(7), mintime(8),     &
                  mintime(6)+mintime(7)+mintime(8)
  WRITE(summary,'("%user: ",f10.2,"  %user: ",f10.2,"   %user: ",f10.2,        &
              &  "   %user: ","       -  ")')                                  &
              (mintime(6)/(mintime(6)+mintime(7)+mintime(8)))*100,             &
              (mintime(7)/(mintime(6)+mintime(7)+mintime(8)))*100,             &
              (mintime(8)/(mintime(6)+mintime(7)+mintime(8)))*100
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              &  "     cpu: ",E10.4)') mintimc(6), mintimc(7), mintimc(8),     &
                  mintimc(6)+mintimc(7)+mintimc(8)
  WRITE(summary,'(" %cpu: ",f10.2,"   %cpu: ",f10.2,                           &
             &"    %cpu: ",f10.2,"    %cpu: ","       -  "/)')                 &
              (mintimc(6)/(mintimc(6)+mintimc(7)+mintimc(8)))*100,             &
              (mintimc(7)/(mintimc(6)+mintimc(7)+mintimc(8)))*100,             &
              (mintimc(8)/(mintimc(6)+mintimc(7)+mintimc(8)))*100

  WRITE(summary,'("COMPARISON OF ELEMENTS_3 _4 AND _5:"/                       &
                 &"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"/)')

  WRITE(summary,'(" Elements_3         Elements_4          Elements_5          &
       &total (c.f. above)"//" MAXIMUMS")')

  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              &  "    user: ",E10.4)') maxtime(4),                             &
                 maxtime(6)+maxtime(7)+maxtime(8), maxtime(9),                 &
                 maxtime(4)+maxtime(6)+maxtime(7)+maxtime(8)+maxtime(9)
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              &  "     cpu: ",E10.4/)') maxtimc(4),                            &
                 maxtimc(6)+maxtimc(7)+maxtimc(8), maxtimc(9),                 &
                 maxtimc(4)+maxtimc(6)+maxtimc(7)+maxtimc(8)+maxtimc(9)

  WRITE(summary,'(" MINIMUMS")')

  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              &  "    user: ",E10.4)') mintime(4),                             &
                 mintime(6)+mintime(7)+mintime(8), mintime(9),                 &
                 mintime(4)+mintime(6)+mintime(7)+mintime(8)+mintime(9)
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              &  "     cpu: ",E10.4/)') mintimc(4),                            &
                 mintimc(6)+mintimc(7)+mintimc(8), mintimc(9),                 &
                 mintimc(4)+mintimc(6)+mintimc(7)+mintimc(8)+mintimc(9)

  WRITE(summary,'("MPI_ALLREDUCE TIME IN BICGSTAB other :"/                    &
                 &"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"/)')

  WRITE(summary,'(" In Ell loop         GCR                 Convergence test   &
       &Total (c.f. above)"//" MAXIMUMS")')

  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              &  "    user: ",E10.4)')                                         &
              maxtime(10)+maxtime(11), maxtime(12),maxtime(13),                &
              maxtime(3)-maxtime(5)-maxtime(6)-maxtime(7)-maxtime(8)-maxtime(9)
  WRITE(summary,'("%user: ",f10.2,"  %user: ",f10.2,"   %user: ",f10.2,        &
              & "   %user: ","       -  ")')                                   &
              ((maxtime(10)+maxtime(11))/(maxtime(3)-maxtime(5)-maxtime(6)     &
              -maxtime(7)-maxtime(8)-maxtime(9)))*100,                         &
              (maxtime(12)/(maxtime(3)-maxtime(5)-maxtime(6)                   &
              -maxtime(7)-maxtime(8)-maxtime(9)))*100,                         &
              (maxtime(13)/(maxtime(3)-maxtime(5)-maxtime(6)-maxtime(7)        &
              -maxtime(8)-maxtime(9)))*100              
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              &  "     cpu: ",E10.4/)')                                        &
              maxtimc(10)+maxtimc(11), maxtimc(12),maxtimc(13),                &
              maxtimc(3)-maxtimc(5)-maxtimc(6)-maxtimc(7)-maxtimc(8)-maxtimc(9)
  WRITE(summary,'(" %cpu: ",f10.2,"  %cpu: ",f10.2,"   %cpu: ",f10.2,          &
              & "   %cpu: ","       -  ")')                                    &
              ((maxtimc(10)+maxtimc(11))/(maxtimc(3)-maxtimc(5)-maxtimc(6)     &
              -maxtimc(7)-maxtimc(8)-maxtimc(9)))*100,                         &
              (maxtimc(12)/(maxtimc(3)-maxtimc(5)-maxtimc(6)                   &
              -maxtimc(7)-maxtimc(8)-maxtimc(9)))*100,                         &
              (maxtimc(13)/(maxtimc(3)-maxtimc(5)-maxtimc(6)-maxtimc(7)        &
              -maxtimc(8)-maxtimc(9)))*100              

  WRITE(summary,'(" MINIMUMS")')

  WRITE(summary,'(" user: ",E10.4,"   user: ",E10.4,"    user: ",E10.4,        &
              &  "    user: ",E10.4)')                                         &
              mintime(10)+mintime(11), mintime(12),mintime(13),                &
              mintime(3)-mintime(5)-mintime(6)-mintime(7)-mintime(8)-mintime(9)
  WRITE(summary,'("%user: ",f10.2,"  %user: ",f10.2,"   %user: ",f10.2,        &
              & "   %user: ","       -  ")')                                   &
              ((mintime(10)+mintime(11))/(mintime(3)-mintime(5)-mintime(6)     &
              -mintime(7)-mintime(8)-mintime(9)))*100,                         &
              (mintime(12)/(mintime(3)-mintime(5)-mintime(6)                   &
              -mintime(7)-mintime(8)-mintime(9)))*100,                         &
              (mintime(13)/(mintime(3)-mintime(5)-mintime(6)-mintime(7)        &
              -mintime(8)-mintime(9)))*100              
  WRITE(summary,'("  cpu: ",E10.4,"    cpu: ",E10.4,"     cpu: ",E10.4,        &
              &  "     cpu: ",E10.4/)')                                        &
              mintimc(10)+mintimc(11), mintimc(12),mintimc(13),                &
              mintimc(3)-mintimc(5)-mintimc(6)-mintimc(7)-mintimc(8)-mintimc(9)
  WRITE(summary,'(" %cpu: ",f10.2,"  %cpu: ",f10.2,"  %cpu: ",f10.2,           &
              & "   %cpu: ","       -  ")')                                    &
              ((mintimc(10)+mintimc(11))/(mintimc(3)-mintimc(5)-mintimc(6)     &
              -mintimc(7)-mintimc(8)-mintimc(9)))*100,                         &
              (mintimc(12)/(mintimc(3)-mintimc(5)-mintimc(6)                   &
              -mintimc(7)-mintimc(8)-mintimc(9)))*100,                         &
              (mintimc(13)/(mintimc(3)-mintimc(5)-mintimc(6)-mintimc(7)        &
              -mintimc(8)-mintimc(9)))*100              

  flops  = REAL(iters*(1+(2*cjiters*ell))*nels_pp*4388.0,iwp)

  WRITE(summary,'(E10.4)') flops

  mflops = flops/(( maxtime(4) + maxtime(6) + maxtime(7) + maxtime(8)   &
                          + maxtime(9) )*1.0E+6 )

  WRITE(summary,'("PE:",I4," Elements_345      ",f6.1," Mflop/s",f6.1," %Peak"/)')&
           numpe, mflops, (mflops/800)*100
 
  flops  = REAL(iters*(cjiters*ell)*nels_pp*4388.0,iwp)
 
  WRITE(summary,'(E10.4)') flops
 
  mflops = flops/(maxtime(7)*1.0E+6)

  WRITE(summary,'("PE:",I4," Matmul Elements_4 ",f6.1," Mflop/s",f6.1," %Peak")')&
           numpe, mflops, (mflops/800)*100
    
END IF

END SUBROUTINE print_timing_info

END MODULE timing
