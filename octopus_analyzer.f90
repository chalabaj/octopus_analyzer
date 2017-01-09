! compile: gfortran -o octopus_analyzer octopus_analyzer_new2.f90
! progam reads line format geometry and transforms it to normal XYZ file 
! energies are extracted and related to the first point to eV units
! neon
!
program octopus_analyze
  implicit none
  INTEGER, external :: Factorial
	INTEGER a,b,i,j,k,l,nlines,nlines_corr,nlines_en,nlines_en_corr,n_folder,step
  INTEGER at,i_atom,Natoms,Ncoors_in_line,atom,Ndists_tot,Ndists_req,i_bond,iost,k_bond
  
	REAL*8, allocatable           :: x(:),y(:),z(:),i_coor(:),i_dist(:),i_dist_req(:),dist(:,:)
	REAL*8, allocatable           :: e_tot(:),ekin_ions(:),e_ion_ion(:),e_electronic(:),e_sum(:)
  REAL*8, allocatable           :: ev_tot(:),evkin_ions(:),ev_ion_ion(:),ev_electronic(:),ev_sum(:)
  INTEGER, allocatable          :: at1(:),at2(:)
  CHARACTER(len=2),allocatable  :: names(:)
  CHARACTER(len=20),allocatable    :: bonds(:)
  
  REAL*8           :: x_dist,y_dist,z_dist
  REAL*8           :: vec1x,vec1y,vec1z,angle,vec2x,vec2y,vec2z,time
	CHARACTER(100)   :: filename,filename3,filename4,outputname,outputname2,outputname3,outputname4,outputname5,arg
	CHARACTER(200)   :: path(20),path2(20),path3(20),path4(20)
	CHARACTER(200)   :: command,command2,command3,command4,command5a,command5b,command6

  CHARACTER(len=5):: str,str2,str3
	REAL*8,parameter :: angsrom=0.52917720859, ev=27.21138602,pi=3.14159265,au_fs=2.418884326505E-2  ! *E-15 = time is in femtoseconds
  LOGICAL(1)          :: f_ex1,f_ex2,f_ex3

!------------------------------------------------------------------------------
  
  n_folder=command_argument_count()
  
  if (n_folder .GT. 20 ) then
    print *,'You can analyze only 20 folder at once'
  STOP 1
  end if
  
  if ( n_folder .EQ. 0) then
      write (*,*)'Please enter a name of a folder with coordinates and energy file:'
      read (*,'(a)')path(1)
      n_folder=1
  else
      j=0  
      do while (j < n_folder)
          j=j+1
          call get_command_argument(j, arg)
	        read(arg,'(A)')path(j)
      end do
  end if 
 
  write(*,*)'How many bonds?'
  read(*, *, IOSTAT=iost)Ndists_req
  if (iost.ne.0)  then 
      write (*,*)'Error: NOT INTEGER VALUE'
      STOP 1
  end if
  
  if (Ndists_req.gt.0)then
      
      allocate( at1(Ndists_req) )
      allocate( at2(Ndists_req) )
      allocate( i_dist_req(Ndists_req) )
      write(*,*)'Please, specify each bond by atom indices, one bond per line.'
  
      do i_bond=1,Ndists_req,1
          read(*,*)at1(i_bond),at2(i_bond)
      end do

  end if
   
  print *,'-----------------------------------------------------------------------------------------------------------------'
!------------------------------------------------------------------------------------FILE HANDLING 
  
  filename=adjustr('/td.general/coordinates')
 	filename3=adjustl(adjustr('/td.general/energy'))
  filename4=adjustl(adjustr('/inp-bohr.xyz'))
 
  do j=1,n_folder,1  

      path(j)=adjustr(path(j))  
	    path2(j)=adjustl(adjustr(path(j)//filename))
      path3(j)=adjustl(adjustr(path(j)//filename3))
      path4(j)=adjustl(adjustr(path(j)//filename4))
    
      outputname=adjustl(path(j)//'-BONDS-MATRIX.dat')
     	outputname2=adjustl(path(j)//'-MOVIE.xyz')
      outputname3=adjustl(path(j)//'-ENERGIES.dat')
      outputname4=adjustl(path(j)//'-ENERGIES_eV.dat')
      outputname5=adjustl(path(j)//'-BONDS.dat')
  
      INQUIRE(FILE = path2(j), EXIST=f_ex1)
      if ( f_ex1 .EQV. .FALSE. ) then
        print *,'Error: coordinates file doesn''t exist in ',trim(path2(j)),'.  SKIPPING'
      	print *,'--------------------------------------------------------------------------------------'
	      GOTO 33
      end if

      INQUIRE(FILE = path3(j), EXIST=f_ex2) 
      if ( f_ex2 .EQV. .FALSE. ) then
        print *,'Error: energy file doesn''t exist in ',trim(path2(j)),'. SKIPPING'
      	print *,'--------------------------------------------------------------------------------------'
	      GOTO 33
      end if

      INQUIRE(FILE = path4(j), EXIST=f_ex3) 
      if ( f_ex3 .EQV. .FALSE. ) then
        print *,'Error: file inp-bohr.xyz doesn''t exist in ',trim(path2(j)),'. SKIPPING'
    	  print *,'--------------------------------------------------------------------------------------'
	      GOTO 33
      end if
     
      print *,f_ex1,f_ex2,f_ex3,'Files check'
      
      print *,'-----------------------------------------------------------------------------------------------------------'

   	  command='wc -l <' // trim(path2(j)) // '> nlines.txt'
      CALL system(command)
      OPEN(101,file='nlines.txt') 
      READ(101,*)nlines 
      close(101)
      nlines_corr=nlines-5

      open(110,file=path4(j),status='OLD')
      read(110,*)Natoms
      read(110,*)
      allocate( names(Natoms) )
      Ncoors_in_line=Natoms*3
      Ndists_tot=Natoms*(Natoms-1)/2
      do i_atom=1,Natoms 
              read(110,*)names(i_atom)
      end do
      close(110)

      allocate( x(Natoms) )
      allocate( y(Natoms) )
      allocate( z(Natoms) )
      allocate( i_coor(Ncoors_in_line))
      allocate( dist(Natoms,Natoms))
      allocate( i_dist(Ndists_tot))
      allocate( bonds(Ndists_tot))
!-----------------------------------------------GEOMETRY CHECK AND TRANSFORMATION	  

      print *,'Reading ',nlines_corr,' geometries from: ',trim(path2(j))
   
      open(102,file=outputname,status='REPLACE')  !BOND-MATRIX  
      open(104,file=outputname2,status='REPLACE') !MOVIE.xyz
      open(111,file=outputname5,status='REPLACE') !requested bonds only  
      open(103,file=path2(j),status='OLD')
!      open(112,file='BONDS-culomn',status='REPLACE')
      read(103,*)
      read(103,*)
      read(103,*)
      read(103,*)
      read(103,*)

!	do i=1,50,1
      do i=1,nlines_corr,1 
        read(103,*)step,time,(i_coor(b),b=1,Ncoors_in_line)
        k_bond=1
        i_atom=1
        do a=0,Ncoors_in_line-3,3
           x(i_atom)=i_coor(a+1)*angsrom
           y(i_atom)=i_coor(a+2)*angsrom
           z(i_atom)=i_coor(a+3)*angsrom
        i_atom=i_atom+1
        end do         
        
        write(104,*)Natoms
        write(104,*)'Step:',Step,'Time (fs):',time   
        
        do l=1,Natoms,1
          write(104,*)names(l),x(l),y(l),z(l)
  	            do k=l+1,Natoms
                    x_dist=(x(l)-x(k))**2
                    y_dist=(y(l)-y(k))**2
                    z_dist=(z(l)-z(k))**2       
                    dist(l,k)=sqrt(x_dist+y_dist+z_dist)
                    i_dist(k_bond)=dist(l,k)
                    
                    if ( i == 1 ) then
                      write(str2,"(I3.2)")l         !,'-',l
                      write(str3,"(I3.2)")k
                      write(str,"(I3.2)") k_bond
                      bonds(k_bond)=trim(adjustl(str))//': '//trim(adjustl(str2))//'-'//trim(adjustl(str3))
                    end if
                    k_bond=k_bond+1
                end do
        end do
        
        time=time*au_fs

        if ( i == 1 ) then
           write(102,*) '#Vazba At1-At2: ',(bonds(k_bond),k_bond=1,Ndists_tot)
           write(102,*) '#Bond are in Angström and and time in femtoseconds'
        end if
        write(102,20)time,(i_dist(k_bond),k_bond=1,Ndists_tot)   
             
        do i_bond=1,Ndists_req,1
            i_dist_req(i_bond)=dist(at1(i_bond),at2(i_bond))   
        end do
        write(111,20) time,(i_dist_req(b),b=1,Ndists_req)
20      format(999F16.8)

      end do
    
    close(102)
	  close(103)
	  close(104)
    close(111)
!ENERGY READING-------------------------------------------------------------------------------------------------------------
	
	  command3='wc -l <' // trim(path3(j)) // '> nlines_en.txt'
	  CALL system(command3)
	  OPEN(105,file='nlines_en.txt') 
	  READ(105,*)nlines_en 
	  close(105)
    nlines_en_corr=nlines_en-5
    
    allocate(e_tot(nlines_corr))
    allocate(ekin_ions(nlines_corr))
    allocate(e_ion_ion(nlines_corr))
    allocate(e_electronic(nlines_corr))
    allocate(e_sum(nlines_corr))
    allocate(ev_tot(nlines_corr))
    allocate(evkin_ions(nlines_corr))
    allocate(ev_ion_ion(nlines_corr))
    allocate(ev_electronic(nlines_corr))
    allocate(ev_sum(nlines_corr))

    print *,'Reading ',nlines_en_corr,' energies from:',path3(j)
  
    open(106,file=path3(j))
  	open(107,file=outputname3,status='REPLACE')
	  open(108,file=outputname4,status='REPLACE')
  
  	write(107,*)'#  t (a.u.)              E_Total (Ha)              E_Kinetic (ions)(Ha)          E_Ion-Ion+Electronic (Ha)'	
  	write(108,*)'#  t (a.u.)             dE_Total (eV)              dE_Kinetic (ions) (eV)        dE_Ion-Ion+Electronic eV)' 	
 
    read(106,*)
    read(106,*)
    read(106,*)
    read(106,*)
    read(106,*)

    do i=1,nlines_en_corr,1  
        if ( i == 1 ) then
           write(107,*) '#Energies are in atomic units and time in femtoseconds.'
           write(108,*) '#Energy differencies are (in electronvols) related to the first step and time is in femtoseconds.'
        end if
        read(106,*)step,time,e_tot(i),ekin_ions(i),e_ion_ion(i),e_electronic(i)
	      time=time*au_fs
        e_sum(i)=e_ion_ion(i)+e_electronic(i)
      
        ev_tot(i)=(e_tot(i)-e_tot(1))*ev
        evkin_ions(i)=(ekin_ions(i)-ekin_ions(1))*ev
        ev_sum(i)=(e_sum(i)-e_sum(1))*ev
	  
        write(107,40)time,e_tot(i),ekin_ions(i),e_sum(i)
 	      write(108,40)time,ev_tot(i),evkin_ions(i),ev_sum(i)
40      format(1F12.6,3E30.14)	
  
   	end do
	
	  close(106)
	  close(107)
 	  close(108)

!END OF ENERGY READING----------------------------------------------------------------------------------------------------

	print *,'Printing results for: ',trim(adjustl(adjustr(path(j))))
	print *,'| All bonds            xmgrace -nxy  ',outputname
	print *,'| Required bonds       xmgrace -nxy  ',outputname5 
	print *,'| Energies             xmgrace -nxy  ',outputname3
	print *,'| dE (in eV)           xmgrace -nxy  ',outputname4
	print *,'| Trajectory           molden        ',outputname2 
	print *,'|------------------------------------------------------------------------------------------------'
  print *,' '

    deallocate( names)
    deallocate( x )
    deallocate( y )
    deallocate( z )
    deallocate( i_coor)
    deallocate( dist)
    deallocate( i_dist)
    deallocate( bonds)
    deallocate( e_tot)
    deallocate( ekin_ions)
    deallocate( e_ion_ion)
    deallocate( e_electronic)
    deallocate( e_sum)
    deallocate( ev_tot)
    deallocate( evkin_ions)
    deallocate( ev_ion_ion)
    deallocate( ev_electronic)
    deallocate( ev_sum)
  
  33 CONTINUE   
  end do
    
end program octopus_analyze

!  INTEGER FUNCTION Factorial(n)
!  IMPLICIT NONE
!  INTEGER, INTENT(IN) :: n
!  INTEGER :: i, Ans 
!  Ans = 1 
!  DO i = 1, n 
!    Ans = Ans * i 
!  END DO 
!  Factorial = Ans 
!  END FUNCTION Factorial
!	command5a=' awk ''{print $2 " " $3 " " $4 " " $5 " " $6}''  '
!	command6=trim(command5a) // ' ' // trim(path3) //  ' > ' // trim(outputname2)
!	CALL system(trim(command6))
!        write(XXX,*)it,(dist(),idist=1,ndist) only required bonds
!    vec1x=x1(i)-x2(i)
!    vec1y=y1(i)-y2(i)
!   vec1z=z1(i)-z2(i)
!    vec2x=x1(i)-x3(i)
!    vec2y=y1(i)-y3(i)
!    vec2z=z1(i)-z3(i)
!  	zlomek=(vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/&
!	  (sqrt(vec1x**2+vec1y**2+vec1z**2)*sqrt(vec2x**2+vec2y**2+vec2z**2))
!  	angle=acos(zlomek)
!    angle=180/pi*angle    
