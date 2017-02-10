 MODULE ANALYSIS
 IMPLICIT NONE
 CONTAINS
 
 SUBROUTINE all_trajs_stat(outputfile,Ifile)
   INTEGER,intent(in)               :: Ifile
   CHARACTER(50),intent(in)         :: outputfile
   INTEGER,parameter                :: Nchannels = 16
   REAL*8,allocatable               :: channel_pop(:,:)

! channel_pop(channel,step) 
 if (Ifile.EQ.1)then

                
 end if
 
 
 END SUBROUTINE all_trajs_stat
 
        SELECT CASE (channel)
         CASE(1)
          channel_pop(channel,step)=channel_pop(channel,step)+1
         CASE(2)
         
         CASE(3)
         
         CASE(4)
         CASE(5)
         CASE(6)
         CASE(7)
         CASE(8)
         CASE(9)
         CASE(10)
         CASE(11)
         CASE(12)
         CASE(13)
         CASE(14)
         CASE(15)  
         CASE(16)       
       END SELECT