#/bin/bash

	
for i in "$@"; do
      input="$i"
      echo "$input"
      nlines=$(wc -l $input | awk {'print $1'})	
Ngeom=1
g=1
      echo $nlines	
      while [ "$Ngeom" -le "$nlines" ];do

        let Ngeom=$g*8000

        if [ "$g" -eq "1" ];then
                echo "Extracting $i in $input,Ngeoms=$Ngeom"
         head -$Ngeom $input | tail -n8 > extracted_$input #new file
        else  
         head -$Ngeom $input | tail -n8 >> extracted_$input
        fi

        let g++
      done

done
