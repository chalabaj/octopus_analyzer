 MODULE ANALYSIS
 IMPLICIT NONE
 CONTAINS
 SUBROUTINE Print_help(err,errfile)
 IMPLICIT  NONE
 integer :: err
 CHARACTER(1000) :: errfile
 
 if (err.EQ.1) then
          write(*,*)'No file to be analyzed or -s -a missing.'
 else if (err.EQ.2)then
    write(*,*)'Couldnt find requested file.',errfile
 end if 
 STOP          
 END SUBROUTINE Print_help
 
 SUBROUTINE channels (dist,Natoms,channel,diss_H,diss_molH)
 IMPLICIT NONE
 REAL*8,  intent(in)           :: dist(:,:)
 INTEGER, intent(in)           :: Natoms 
 INTEGER, intent(out)          :: channel,diss_H,diss_molH ! output vars
 INTEGER                       :: hh
 INTEGER, allocatable          :: diss_H_at(:)
 REAL*8, parameter             :: H_diss_limit = 3.000, HH_bond_dist = 1.5
  
! WAT DIM ANALYSIS at each timestep
! 1: nothing, water monomer
! 2: 1H diss 
! 3: 2xH_at
! 4: H2 mol

        diss_H = 0
        diss_molH = 0
        channel = 1
        allocate ( diss_H_at(4) )
        
        do hh=2,Natoms,1   ! 1 is oxygen atom
          if ( dist(1,hh).gt.H_diss_limit.AND.dist(2,3).gt.HH_bond_dist ) then
            diss_H = diss_H + 1
            diss_H_at(dissH) = hh          ! which H(hh) is dissociated
          ! if not dissociated then where the hydrogen is? O1 or O2 1
          end if
        end do

        if ( diss_H.EQ.0 ) then 
          channel = 1   !         1: Nothing happens, water monomer     
        else if ( diss_H.EQ.1 ) then 
          channel = 2   !         2: 1H diss                             
        else if ( diss_H.EQ.2 ) then         
                if ( dist(2,3).gt.HH_bond_dist ) then
                 channel = 3 !    3:  2xH_at
                else if ( dist(2,3).le.HH_bond_dist ) then
                 diss_molH = 1
                 channel = 4 !    4:  H2 mol
                end if
        end if
        
        deallocate ( diss_H_at ) ! next geometry starts from begining
        
!-------------------cluster analysis done ------------------------------------------  
 
 END SUBROUTINE channels 
 END MODULE ANALYSIS
