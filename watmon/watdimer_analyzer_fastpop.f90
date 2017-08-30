! gfortran -o watdim_analyzer  analysis.f90 watdimer_analyzer_fastpop.f90  -g -fcheck=all -Wall

 program  water_dimer_analysis
 USE ANALYSIS
 IMPLICIT NONE
 
 REAL*8,allocatable       :: channel_pop(:,:),times(:),totpop(:)
 REAL*8, allocatable      :: x(:),y(:),z(:),dist(:,:)  
 REAL*8                   :: x_dist,y_dist,z_dist
 REAL*8                   :: time,OOdist,tim,dt
 REAL*8,parameter         :: au_fs=2.418884326505E-2
 CHARACTER(3)             :: mode
 CHARACTER(50)            :: inputfile(1000),arg
 CHARACTER(50)            :: outputfile
 CHARACTER(100)           :: tt,command
 CHARACTER, allocatable   :: names(:)
 
 INTEGER                  :: dissH,dissmolH,ho1,ho2  ! water dimer analysis variables 
 INTEGER                  :: outlines,istep,maxistep !maxistep is max step in output file
 INTEGER                  :: Nchannels, Maxstep !simulation Maxsteps
 INTEGER                  :: igeom,Ngeoms,Natoms,i_atom,nlines,channel
 INTEGER                  :: l,k,j,step,Nargs  ! loops
 
 LOGICAL                  :: f_ex1
!##############################################################################
 Nchannels = 16
 Maxstep = 400000
 dt = 4 * au_fs! timestpe for SH
 outputfile = 'dataall.dat'
 
 ! Error handling--------------------------------------------------------------
 Nargs=command_argument_count()
 call get_command_argument(1, mode)
  if (mode.EQ.'-st')then
 GOTO 12
 else if ( Nargs.LT.2 )then
     call Print_help(1,'')
 end if


!Reading input files-----------------------------------------------------------
 j=2               ! j-th geometry
      do while (j.LE.Nargs)
       call get_command_argument(j, arg)  ! first argument is type of movie: octopus or SH
       read(arg,'(A)')inputfile(j)       
       INQUIRE(FILE =inputfile(j), EXIST=f_ex1)
       if ( f_ex1 .EQV. .FALSE. ) then
          call Print_help(2,inputfile(j))
       end if
       j=j+1
      end do

!Geometry###################################################################### 
!----------proces movie------------------------- 
 
   print *,'starting analysis for ',Nargs-1,' files.'
     
   open(111,file=outputfile,status='REPLACE',access='append') 
   write(111,*) '#Time,     Step,  Channel,   O-O dist,     xH-O1, xH-O1, diss H, diss mol H2'
   close(111)
   
   do j=2,Nargs,1
   print *,'Analysis file: ', inputfile(j)
        command='wc -l <' // inputfile(j) // '> nlines.txt'
        CALL system(command)
        OPEN(101,file='nlines.txt') 
        READ(101,*)nlines 

        close(101)
 !       CALL system('rm nlines.txt')
    
        open(110,file=inputfile(j),status='OLD')
          

        open(111,file=outputfile,status='OLD',access='append') 
        read(110,*)Natoms
        REWIND(110) 
        allocate( x(Natoms) )
        allocate( y(Natoms) )
        allocate( z(Natoms) )
        allocate( names(Natoms) )
        allocate( dist(Natoms,Natoms))

        Ngeoms = nlines / (Natoms+2)  
        print *,'Nlines, Natoms, Ngeoms',Nlines,'  ',Natoms,'  ',Ngeoms
        print *,'--------------------------------------------'
        do igeom=1,Ngeoms,1


           read(110,*)Natoms
           if (mode.EQ.'-oc')then
             read(110,*)tt,step,tt,tt,time  !for testing geometrie added tt
           else if (mode.EQ.'-sh')then
             read(110,*)tt,tt,step
             time = step * dt ! dt from input.in
             step = step - 1  ! octopus starts from 0 but SH numbering from 1
           end if
           
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

            write(111,5)time,step,channel,dist(1,2),ho1,ho2,dissH,dissmolH
5           format(1F9.4,2I8.1,1F16.8,4I8.1) 
        end do  
              
      close(110)
      close(111)  
        deallocate( x )
        deallocate( y )
        deallocate( z )
        deallocate( names )
        deallocate( dist )   
  end do

!-----calling stastics----------------

12   print *,'processing populations'
   
   command='wc -l <' // outputfile // '> outlines.txt'
    CALL system(command)
    OPEN(113,file='outlines.txt') 
    READ(113,*)outlines 
    close(113)

! array allocation and filling it with zero population    
    allocate ( channel_pop(1:Nchannels,0:Maxstep) ) 
    allocate ( times(0:Maxstep) )   
    allocate ( totpop(0:Maxstep))
!    channel_pop(0,0) = 0
!    totpop(0) = 0
    do channel=1,Nchannels,1
     do step=0,Maxstep,1
       channel_pop(Nchannels,step) = 0
       totpop(step) = 0
     end do
    end do
    
    step = 0
    channel = 1
    maxistep = 1
!-----------READING DATAALL    
    open(112,file=outputfile,status='OLD')
    read(112,*) !first line with comments 
    do istep=1,outlines-1,1  !through entire file except first line
        read(112,*)tim,step,channel,OOdist,ho1,ho2,dissH,dissmolH
        channel_pop(channel,step)=channel_pop(channel,step)+1 
        times(step) = tim
        if (step.GT.maxistep)then !setting maximum number of step for sumilation
          maxistep = step 
        end if 
    end do
    close(112)
    
!-----------POPULATIONS      
    open(114,file='RESULTSTRAJS.dat',status='REPLACE')        
    write(114,*)'#times,  totpop,  channels....'
     do step=0,maxistep,1  !only through maxistep since simualations possible didnt finished
      
       do channel=1,Nchannels,1
        totpop(step)= totpop(step)+channel_pop(channel,step)
       end do 
       
       ! after sum need normed population, cant put upward since i cant do the norm pop before sum
       do channel=1,Nchannels,1 
        channel_pop(channel,step) = channel_pop(channel,step)/totpop(step)*100 ! print %
       end do
!       write(114,6)times(step),totpop(step),(channel_pop(channel,step),channel=1,Nchannels)
            
       write(114,6)times(step),(channel_pop(channel,step),channel=1,Nchannels)
6      format(18F16.8)

     end do
     close(114)
 
    deallocate ( channel_pop ) 
    deallocate ( times )
    STOP

 end program water_dimer_analysis
 
 