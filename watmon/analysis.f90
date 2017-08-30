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
 
 SUBROUTINE channels (dist,Natoms,channel,dissH,dissmolH,ho1,ho2)
 IMPLICIT NONE
 REAL*8,  intent(in)           :: dist(:,:)
 INTEGER, intent(in)           :: Natoms 
 INTEGER, intent(out)          :: channel,dissH,dissmolH,ho1,ho2  ! output vars
 INTEGER                       :: Ndiss_H,h1,h2,diss_molH,hh,h_o1,h_o2
 INTEGER, allocatable          :: diss_H_at(:)
 REAL*8, parameter             :: OH_bond_dist = 2.000 ,OO_bond_dist = 4.500, H_diss_limit = 3.000
 REAL*8, parameter             :: HH_bond_dist = 1.5
  
! WAT DIM ANALYSIS at each timestep
! diss_H - dissociated H atoms, if diss_H > 2 check for molecular hydrogen 
! h_o1 h_o2 number of hydrogen atoms on each oxygen   
! channels: 
! 1 H2O...H2O bonded
! 2 H2O+H2O disscociated
! 3 PT H3O...OH bonded
! 4 PT H3O+OH  dissociated
! 5 H diss H20...OH bonded
! 6 H diss H20+OH dissociated

! 7  2x at H diss,  O + H2O dissociated 
! 8  2x at H diss,  O...H2O bonded  
! 9  2x at H2 diss, 2x OH dissociated  
! 10 2x at H2 diss, OH...OH bondedg

! 11 1 mol H2 diss, O + H2O dissociated 
! 12 1 mol H2 diss, O...H2O bonded 
! 13 1 mol H2 diss, 2x OH dissociated  
! 14 1 mol H2 diss, OH...OH bondedg

        h_o1 = 0
        h_o2 = 0
        Ndiss_H = 0
        diss_molH = 0
        channel = 1
        allocate ( diss_H_at(4) )
        
        do hh=3,Natoms,1   ! 1 and 2 are oxygen atoms
          if ( dist(1,hh).gt.H_diss_limit.AND.dist(2,hh).gt.H_diss_limit ) then
            Ndiss_H = Ndiss_H + 1
            diss_H_at(Ndiss_H) = hh          ! which H(hh) is dissociated
          ! if not dissociated then where the hydrogen is? O1 or O2 
          else if( dist(1,hh).lt.dist(2,hh).AND.dist(1,hh).lt.OH_bond_dist )then 
          h_o1 = h_o1 + 1
          else if( dist(2,hh).lt.dist(1,hh).AND.dist(2,hh).lt.OH_bond_dist ) then
          h_o2 = h_o2 + 1
          end if
        end do

        if ( Ndiss_H.EQ.0 ) then 
! NO PT or H diss  
           if ( h_o1.eq.2.AND.h_o2.eq.2 ) then
                if ( dist(1,2).gt.OO_bond_dist ) then
                  channel = 2 ! 2 H2O+H2O disscociated 
                else 
                 channel = 1 ! 1 H2O...H2O bonded 
                end if
           else if ( h_o1.eq.3.OR.h_o2.eq.3 ) then
                if ( dist(1,2).gt.OO_bond_dist ) then
                  channel = 4   ! 4 PT H3O+OH  dissociated
                else
                  channel = 3   ! 3 PT H3O...OH bonded  
                end if 
           end if 
                  
        else if ( Ndiss_H.EQ.1 ) then 

                if ( dist(1,2).gt.OO_bond_dist ) then
                  channel = 6   ! 6 H diss, H20+OH dissociated
                else
                  channel = 5   ! 5 H diss, H20...OH bonded 
                end if 
                
                if ( h_o1.eq.3.OR.h_o2.eq.3 ) then
                if ( dist(1,2).gt.OO_bond_dist ) then
                  channel = 15   ! 15 PT H3O+O + H  dissociated TEST
                else
                  channel = 16   ! 16 PT H3O...O + H bonded  TEST
                end if 
           end if 
                               
        else if ( Ndiss_H.GE.2 ) then         
                do h1=1,Ndiss_H,1
                    do h2=h1+1,Ndiss_H,1
                       if ( dist(diss_H_at(h1),diss_H_at(h2)).le.HH_bond_dist ) then
                          diss_molH = diss_molH + 1
                       end if
                    end do
                 end do
                 if ( diss_molH.eq.0 )then ! not mol H2 but 2x H
                     if ( h_o1.eq.0.OR.h_o2.eq.0 ) then
                         if ( dist(1,2).gt.OO_bond_dist ) then
                            channel = 7     ! 7 2x at H diss, O + H2O dissociated 
                         else   
                            channel = 8     ! 8 2x at H diss, O...H2O bonded 
                         end if
                     else if ( h_o1.eq.1.AND.h_o2.eq.1 ) then
                         if ( dist(1,2).gt.OO_bond_dist ) then
                            channel = 9     ! 9 2x at H2 diss, 2x OH dissociated  
                         else
                            channel = 10    ! 10 2x at H2 diss, OH...OH bonded
                         end if 
                     end if   
                 else if ( diss_molH.eq.1 )then
                     if ( h_o1.eq.0.OR.h_o2.eq.0 ) then
                         if ( dist(1,2).gt.OO_bond_dist ) then
                            channel = 11     ! 11 1 mol H2 diss, O + H2O dissociated 
                         else   
                            channel = 12     ! 12 1 mol H2 diss, O...H2O bonded 
                         end if
                     else if ( h_o1.eq.1.AND.h_o2.eq.1 ) then
                         if ( dist(1,2).gt.OO_bond_dist ) then
                            channel = 13     ! 13 mol H2 diss, 2x OH dissociated  
                         else
                            channel = 14    ! 14 mol H2 diss, OH...OH bonded
                         end if 
                     end if      
                 end if
                 
        else if ( Ndiss_H.Gt.2 ) then
         channel = 10000  ! too much to think about all the channels
        end if
        dissH = Ndiss_H
        dissmolH = diss_molH
        ho1 = h_o1
        ho2 = h_o2

!        print *, channel,h_o1,h_o2,dist(1,2),Ndiss_H,diss_molH

        deallocate ( diss_H_at ) ! next geometry starts from begining
        
!-------------------cluster analysis done ------------------------------------------  
 
 END SUBROUTINE channels 
 END MODULE ANALYSIS
