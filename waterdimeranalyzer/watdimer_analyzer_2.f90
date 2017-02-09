 program  water_dimer_analysis
 USE ANALYSIS
 IMPLICIT NONE

 REAL*8, allocatable          :: x(:),y(:),z(:),dist(:,:)  
 REAL*8                       :: x_dist,y_dist,z_dist,time
 INTEGER                      :: igeom,Ngeoms,Natoms,i_atom,Nfiles,nlines,channel
 INTEGER                      :: l,k,j  ! loops
 CHARACTER(50)                :: inputfile(1000),arg,outputfile(1000)
 CHARACTER(100)               :: tt,path(1000),command
 CHARACTER, allocatable       :: names(:)
 LOGICAL                      :: f_ex1
 
 INTEGER                      :: dissH,dissmolH,ho1,ho2                ! water dimer analysis variables


 ! Error handling-----------------------------------------------------------------
 Nfiles=command_argument_count()
 if ( Nfiles.LT.1 )then
  call Print_help(1,'')
 end if 

 !Reading input files-------------------------------------------------------------
 
 j=1               ! j-th geometry
      do while (j.LE.Nfiles)
          call get_command_argument(j, arg)     ! j+1 because first argument is controling the units and the rest are the files to process  
	        read(arg,'(A)')inputfile(j)       
          INQUIRE(FILE =inputfile(j), EXIST=f_ex1)
          if ( f_ex1 .EQV. .FALSE. ) then
            call Print_help(2,path(j))
          end if
          outputfile(j)='RESULTS-'//inputfile(j)
          j=j+1
      end do
! channel = 0  ! noreaction at the beginning
 
 
!Geometry-------------------------------------------------------------
  do j=1,Nfiles,1
     
        command='wc -l <' // inputfile(j) // '> nlines.txt'
        CALL system(command)
        OPEN(101,file='nlines.txt') 
        READ(101,*)nlines 
        close(101)
 !       CALL system('rm nlines.txt')
    
        open(110,file=inputfile(j),status='OLD')  
        open(111,file=outputfile(j),status='REPLACE') 
        write(111,*) '#Time,         Channel,  O-O dist,       xH-O1, xH-O1, diss H, diss mol H2'
        read(110,*)Natoms
        REWIND(110) 
        allocate( x(Natoms) )
        allocate( y(Natoms) )
        allocate( z(Natoms) )
        allocate( names(Natoms) )
        allocate( dist(Natoms,Natoms))

        Ngeoms = nlines / (Natoms+2)  
      
        do igeom=1,Ngeoms,1


           read(110,*)Natoms
           read(110,*)tt,tt,tt,tt,tt,time  !for testing geometrie added tt, removie for movie so it is possible to read time
           do i_atom=1,Natoms,1
               read(110,*)names(i_atom),x(i_atom),y(i_atom),z(i_atom) 
           end do
      
           do l=1,Natoms,1
                do k=l+1,Natoms,1               ! distance matrix
                  x_dist=(x(l)-x(k))**2
                  y_dist=(y(l)-y(k))**2
                  z_dist=(z(l)-z(k))**2       
                  dist(l,k)=sqrt(x_dist+y_dist+z_dist)
                end do
           end do
           
           call channels(dist,Natoms,channel,dissH,dissmolH,ho1,ho2)

            write(111,5)time,channel,dist(1,2),ho1,ho2,dissH,dissmolH
5           format(1F12.8,1I8.1,1F16.8,4I8.1) 
        end do  
              
      close(110)
      close(111)  
        deallocate( x )
        deallocate( y )
        deallocate( z )
        deallocate( names )
        deallocate( dist )   
  end do
  
 end program water_dimer_analysis
 
 