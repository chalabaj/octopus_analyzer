 program  water_dimer_analysis
 USE ANALYSIS
! USE STATISTICS
 IMPLICIT NONE

 REAL*8, allocatable      :: x(:),y(:),z(:),dist(:,:)  
 REAL*8                   :: x_dist,y_dist,z_dist,time,OOdist,tim
 INTEGER                  :: igeom,Ngeoms,Natoms,i_atom,Nfiles,nlines,channel
 INTEGER                  :: l,k,j,step,Nargs  ! loops
 CHARACTER(50)            :: inputfile(1000),arg
 CHARACTER=               :: outputfile = 'data_all.dat'
 CHARACTER(100)           :: tt,path(1000),command
 CHARACTER, allocatable   :: names(:)
 
 INTEGER                  :: outlines,istep,maxistep !maxistep is max step in output file
 INTEGER,parameter        :: Nchannels = 16, Maxstep = 400000
 REAL*8,allocatable       :: channel_pop(:,:),times()
 LOGICAL                  :: f_ex1
 CHARACTER(2)             :: mode
 
 INTEGER                  :: dissH,dissmolH,ho1,ho2                ! water dimer analysis variables


 ! Error handling-----------------------------------------------------------------
 Nargs=command_argument_count()
 if ( Nfiles.LT.2 )then
  call Print_help(1,'')
 end if 

 !Reading input files-------------------------------------------------------------
 
 j=1               ! j-th geometry
      do while (j.LE.Nargs)
          call get_command_argument(j+1, arg)  ! first argument is operation mode
	        read(arg,'(A)')inputfile(j)       
          INQUIRE(FILE =inputfile(j), EXIST=f_ex1)
          if ( f_ex1 .EQV. .FALSE. ) then
            call Print_help(2,path(j))
          end if
          j=j+1
      end do
   
! channel = 0  ! noreaction at the beginning

!-----calling stastics----------------
 Nfiles=command_argument_count()-1 
 
 call get_command_argument(1, mode)
  if (mode.EQ.'-s') then
    command='wc -l <' // outputfile // '> outlines.txt'
    CALL system(command)
    OPEN(113,file='outlines.txt') 
    READ(113,*)outlines 
    close(113)
    CALL system('rm outlines.txt')

! array allocation and filling it with zero population    
    allocate ( channel_pop(Nchannels,Maxstep) ) 
    allocate ( times(Maxstep) )   
    
    do channel=1,Nchannels,1
     do step=1,Maxstep
       channel_pop(Nchannels,Maxstep) = 0
     end do
    end do
    
    open(112,file=outputfile,status='OLD') 
     
     totpop(step) = 0
     step = 0
     channel = 1
     maxistep = 1
     do istep=1,outlines,1  !through entire file
        read(112,5)tim,step,channel,OOdist,ho1,ho2,dissH,dissmolH 
        channel_pop(channel,step)=channel_pop(channel,step)+1 
        times(step) = tim
        if (step.GT.maxistep)then
          maxistep = step 
        end if 
     end do
     
     close(112)
       
     do step=1,maxistep,1  !only through maxistep since simualations possible didnt finished
      
       do channel=1,Nchannels,1
        totpop(step)= totpop(step)+channel_pop(channel,step))
       end do 
       
       ! after sum need normed population, cant put upward since i cant do the norm pop before sum
       do channel=1,Nchannels,1 
        channel_pop(channel,step) = channel_pop(channel,step)/totpop(step)
       end do
     
       write(114,6)times(step),(channel_pop(channel,step),channel=1,Nchannels)
6      format(17F16.8)

     end do

 
  else if (mode.EQ.'-a')GOTO 11
  end if
  
!Geometry------------------------------------------------------------- 
!----------proces movie-------------------------  
11  do j=1,Nfiles,1
  
        command='wc -l <' // inputfile(j) // '> nlines.txt'
        CALL system(command)
        OPEN(101,file='nlines.txt') 
        READ(101,*)nlines 
        close(101)
 !       CALL system('rm nlines.txt')
    
        open(110,file=inputfile(j),status='OLD')  
        open(111,file=outputfile,status='OLD',access='append') 
        write(111,*) '#Time,   Step,      Channel,  O-O dist,       xH-O1, xH-O1, diss H, diss mol H2'
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
           read(110,*)tt,step,tt,tt,time  !for testing geometrie added tt, removie for movie so it is possible to read time
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

            write(111,5)time,step channel,dist(1,2),ho1,ho2,dissH,dissmolH
5           format(1F12.8,2I8.1,1F16.8,4I8.1) 
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
 
 