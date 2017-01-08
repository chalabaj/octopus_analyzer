!gfortran -o hbond_analyzer hbondanalyzer.f90

	program readall
	integer j,i,natom,step,nlines,nlines_corr
	character atom(1)
	real x1(1000000),y1(1000000),z1(100000),x2(1000000),y2(1000000),z2(100000),x3(1000000),y3(1000000),z3(100000),angsrom
	real distx12,disty12,distz12,distx13,disty13,distz13,dist12,dist13,junk
        real*8  :: vec1x,vec1y,vec1z,angle,pi=3.14159265
        real*8  :: vec2x,vec2y,vec2z,cit,jmen,zlomek
	character(8)   :: junkspace,resultname
	CHARACTER(100) :: filename,filename2,outputname,outputname2,outputname3
	CHARACTER(200) :: path,path2,path3
	CHARACTER(200) :: command,command2,command3,command4,command5a,command5b,command6

	angsrom=0.52917720859

!------------------------------------------------------------------------------
	write (*,*)'Please enter a name of a folder with coordinate and energy file:'
	read (*,'(a)')path

	path=trim(adjustr(adjustl(path)))
	print *,'-------------------------------------------------------'

	filename=adjustl(adjustr('/td.general/coordinates'))
	
	path2=adjustl(adjustr(path//filename))
	print *,'Reading file:',path2

	command='wc -l <' // trim(path2) // '> nlines.txt'

!	print *,command

	CALL system(command)

	outputname=adjustl(path//'-RESULTS.dat')
	outputname3=adjustl(path//'-MOVIE.xyz')

!-------------------------------------------------------------------------------------	
	OPEN(101,file='nlines.txt') 
	READ(101,*)nlines 
	close(101)
	nlines_corr=nlines-5
!	print *,'Number of geometries:'
!	print *,nlines_corr
	
        open(102,file=outputname,status='REPLACE')
        open(103,file=path2)
	open(104,file=outputname3)
	read(103,*)
        read(103,*)
        read(103,*)
        read(103,*)
        read(103,*)

	do i=1,nlines,1 
 
        read(103,*,end=10)step,junk,x1(i),y1(i),z1(i),x2(i),y2(i),z2(i),x3(i),y3(i),z3(i)

	x1(i)=x1(i)*angsrom
        y1(i)=y1(i)*angsrom
        z1(i)=z1(i)*angsrom
        x2(i)=x2(i)*angsrom
        y2(i)=y2(i)*angsrom
        z2(i)=z2(i)*angsrom
        x3(i)=x3(i)*angsrom
        y3(i)=y3(i)*angsrom
        z3(i)=z3(i)*angsrom
	

        write(104,*)"3"
        write(104,*)i
        write(104,*)"O",x1(i),y1(i),z1(i)
        write(104,*)"H",x2(i),y2(i),z2(i)
        write(104,*)"H",x3(i),y3(i),z3(i)



	distx12=(x1(i)-x2(i))**2
        disty12=(y1(i)-y2(i))**2
        distz12=(z1(i)-z2(i))**2
        
        dist12=sqrt(distx12+disty12+distz12)
!	print *,'OH1',dist12

	distx13=(x1(i)-x3(i))**2
        disty13=(y1(i)-y3(i))**2
        distz13=(z1(i)-z3(i))**2
        dist13=sqrt(distx13+disty13+distz13)
!	print *,'OH2',dist13

 
        vec1x=x1(i)-x2(i)
        vec1y=y1(i)-y2(i)
        vec1z=z1(i)-z2(i)
        vec2x=x1(i)-x3(i)
        vec2y=y1(i)-y3(i)
        vec2z=z1(i)-z3(i)

	zlomek=(vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/&
	(sqrt(vec1x**2+vec1y**2+vec1z**2)*sqrt(vec2x**2+vec2y**2+vec2z**2))
	angle=acos(zlomek)
        angle=180/pi*angle
	write(102,20)step,dist12,dist13,angle
20      format(I6,F12.6,F12.6,F12.5)

	end do


10	continue
	close(102)
	close(103)
	close(104)
!-----------------------------------------------------------------
	outputname2=adjustl(path//'-ENERGIES.dat')

	filename2=adjustl(adjustr('/td.general/energy'))
	path3=adjustl(adjustr(path//filename2))

	command5a=' awk ''{print $1 " " $3 " " $4 " " $6}''  '
	command6=trim(command5a) // ' ' // trim(path3) //  ' > ' // trim(outputname2)
	print *,command6

	CALL system(trim(command6))
	print *,'---------------------------------------------------------------------'
	print *,'| Geometrical parameters were saved to:   xmgrace -nxy  ',outputname
	print *,'| Energies were saved to:                 xmgrace -nxy  ',outputname2
	print *,'| Trajectory was  saved to:               molden ',outputname3 
	print *,'---------------------------------------------------------------------'


!--------------------------------------------------------------------
	end program readall
