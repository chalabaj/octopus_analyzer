!111
! compile: gfortran -o octopus_analyzer octopus_analyzer.f90
! progam reads line format geometry and transforms it to normal XYZ file 
! energies are extracted and related to the first point to eV units
program octopus_analyze
  implicit none
	INTEGER a,b,i,j,k,l,nlines,nlines_corr,nlines_en,nlines_en_corr,n_folder,step
  INTEGER at,i_atom,Natoms,Ncoors_in_line,atom,Ndists_tot,Ndists_req,i_bond,iost,k_bond
  
  INTEGER hh,h_o1,h_o2,Ndiss_H,diss_H_at(:),h1,h2,channel                              ! water dimer analysis variables
  
	REAL*8, allocatable           :: x(:),y(:),z(:),i_coor(:),i_dist(:),i_dist_req(:),dist(:,:)         ! bond analysis variables
	REAL*8, allocatable           :: e_tot(:),ekin_ions(:),e_ion_ion(:),e_electronic(:),e_sum(:)        ! energy analysis variables
  REAL*8, allocatable           :: ev_tot(:),evkin_ions(:),ev_ion_ion(:),ev_electronic(:),ev_sum(:)   ! energy analysis variables
  INTEGER, allocatable          :: at1(:),at2(:)                                                      ! atom's label for selected bonds analysis
  CHARACTER(len=2),allocatable  :: names(:)                                                           ! atom's names
  CHARACTER(len=20),allocatable :: bonds(:)
  
  REAL*8           :: x_dist,y_dist,z_dist,time
	CHARACTER(100)   :: filename,filename3,filename4,outputname,outputname2,outputname3,outputname4,outputname5,arg
	CHARACTER(200)   :: path(1000),path2(1000),path3(1000),path4(1000)
	CHARACTER(200)   :: command,command2,command3,command4,command5a,command5b,command6                  ! calling linux comand

  CHARACTER(len=5) :: str,str2,str3                                                                    ! text explaining bonds in BOND-MATRIX
	REAL*8,parameter :: angsrom=0.52917720859, ev=27.21138602,pi=3.14159265,au_fs=2.418884326505E-2      ! *E-15 = time is in femtoseconds
  LOGICAL(1)       :: f_ex1,f_ex2,f_ex3                                                                ! file exists status

!------------------------------------------------------------------------------
  
  n_folder=command_argument_count()
  
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

! READING NAME OF ATOMS FROM 
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
      open(103,file=path2(j),status='OLD')        !
      read(103,*)
      read(103,*)
      read(103,*)
      read(103,*)
      read(103,*)

      do i=1,nlines_corr,1                                               !main loop for reading and ANALYSING GEOMETRY
        read(103,*)step,time,(i_coor(b),b=1,Ncoors_in_line)                  
        time=time*au_fs
             
!  COORDINATE TRANSFORM {line to collumn}
        i_atom=1
        do a=0,Ncoors_in_line-3,3
           x(i_atom)=i_coor(a+1)*angsrom
           y(i_atom)=i_coor(a+2)*angsrom
           z(i_atom)=i_coor(a+3)*angsrom
        i_atom=i_atom+1
        end do         
 
! MOVIE       
       write(104,*)Natoms
       write(104,*)'Step:',Step,'Time (fs):',time                       !heading for each step in movie
        
        k_bond=1
        do l=1,Natoms,1
            write(104,*)names(l),x(l),y(l),z(l) ! writing XYZ to movie file
            do k=l+1,Natoms               ! distance matrix
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

! WAT DIM ANALYSIS at each timestep
! diss_H - dissociated H atoms, if diss_H > 2 check for molecular hydrogen 
! h_o1 h_o2 number of hydrogen atoms on each oxygen   
! channels: 
! 1 H2O...H2O bonded
! 2 H2O+H2O disscociated
! 3 PT H3O...OH bonded
! 3 PT H3O+OH  dissociated
! 5 H diss H20 OH
        h_o1 = 0
        h_o2 = 0
        Ndiss_H = 0
        diss_molH = 0
        channel = 1
        allocate ( diss_H_at(4) )
        
        do hh=3,Natoms,1
          if ( dist(1,hh).gt.3.AND.dist(2,hh).gt.3 ) then
            Ndiss_H = Ndiss_H + 1
            diss_H_at(Ndiss_H) = hh          ! which H(hh) is dissociated
          else if ( dist(1,hh).lt.dist(2,hh).AND.dist(1,hh).lt.1.4 ) then   ! if not dissociated then where the hydrogen is? O1 or O2
 	           h_o1=h_o1+1
          else if ( dist(2,hh).lt.dist(1,hh).AND.dist(2,hh).lt.1.4 ) then
             h_o2=h_o2+1
          end if
        end do

          
        if ( Ndiss_H.EQ.0 ) then 
! 2 H2O+H2O disscociated  // 1 H2O...H2O bonded     
          if ( h_o1.eq.2.AND.h_o2.eq.2 ) then
                if ( dist(1,2).gt.4.0 ) then
                  channel = 2
                else
                  channel = 1 
                end if
          GOTO 5
! 3 PT H3O...OH bonded   
          else if ( h_o1.eq.3.OR.h_o2.eq.3 ) then
                if ( dist(1,2).gt.4.0 ) then
                  channel = 3   ! 3 PT H3O+OH  dissociated
                else
                  channel = 4
                end if 
          end if 
                  
        else if ( Ndiss_H.EQ.1 ) then 
        XXX
        else if ( Ndiss_H.GE.2 ) then         
          do h1=1,Ndiss_H,1
            do h2=h1+1,Ndiss_H,1
               if ( dist(diss_H_at(h1),diss_H_at(h2)).le.1.5 ) then
                diss_molH = diss_molH + 1
               end if
             end do
           end do
        end if


        end if

! molecular hydrogen        

5         
        deallocate ( diss_H_at ) ! next geometry starts from begining
        
!-------------------cluster analysis done ------------------------------------------        
        if ( i == 1 ) then
           write(102,*) '#Bond are in Angstroem and time in femtoseconds'
           write(102,*) '#Time i-th bond: Atom1-Atom2: ',(bonds(k_bond),k_bond=1,Ndists_tot)
        end if
        write(102,20)time,(i_dist(k_bond),k_bond=1,Ndists_tot)   
             
        do i_bond=1,Ndists_req,1
            i_dist_req(i_bond)=dist(at1(i_bond),at2(i_bond))   
        end do
        write(111,20)time,(i_dist_req(b),b=1,Ndists_req)
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

