#/bin/bash

file=$1

# Combines the samples from the MCMC from a single latent 
# sample from each gene. This reduces the size of the files
# that need to be read into R.


# Number of MCMC files
numFiles=$(head -1 ${file}| sed 's/,/ /g'| wc -w)

# Number of latent sequences plus the header
numLines=$(wc -l $file | awk '{print $1}')

for ((i=1;i<=numFiles;i++))
do
	name=$(head -1 ${file}|awk -F "," -v col=$i '{printf $col}')
	outfile=${name}/output

	grep "Date range" ${outfile} | awk '{print $3 $10}' >> dateRange
done

tipDateCheck=$(sort dateRange | uniq |wc| awk '{print $1}')

if [[ $tipDateCheck != 1 ]] 
then

	echo All runs of HIVtree must be run with the same tipDate time unit and have the same last sample date. 
	echo This information can be checked on the line starting with "Date range" in the output to screen of the MCMC. 
	echo See the manual for more information
else 

	time=$(head -1 dateRange | awk -F , '{print $2}')
	lastSample=$(head -1 dateRange | awk -F , '{print $1}' | sed 's/(//g')
	echo The tipDate time unit is ${time} 
	echo The time of the last sample is ${lastSample}.
fi

rm dateRange

# For each row in the csv file
for ((line=2;line<=numLines;line++));
do

	echo Processing row $line in the csv file.
	# Name the output file based on the first sequence name in a row

	#seq1=$(tail -n+$line $file |head -1 |awk -F "," -v col=1 '{printf $col}')
	lineText=$(tail -n+$line $file |head -1)

	# The first sequence is not missing
	pat='^([^,]+),(.*\b)'

	# The first sequence is missing
	pat1='(,+)([^,]+)(,?)(.*\b)'

	# Finds the name of the first sequence in a row
	if [[ "$lineText" =~ $pat ]]; then
    		seq1=${BASH_REMATCH[1]}

	elif [[ "$lineText" =~ $pat1 ]]; then
    		seq1=${BASH_REMATCH[2]}
	else

		echo Check file format. The file should be in csv format. Problem in row 
		echo $lineText
		echo Exiting. 
		exit
	fi 

	touch ${seq1}.txt

	# For each column in the row
	for ((i=1;i<=numFiles;i++));
	do
		
		# Find the sequence name
		seq=$(tail -n+$line $file |head -1 |awk -F "," -v col=$i '{printf $col}')

		# This latent sequence does not exist for this gene
		if [ -z "$seq" ] 
		then

			# Name the column in the output the MCMC file that does not have
			# the latent sequence. This will be a blank column except the header
			head -1 $file | awk -F "," -v col=$i '{printf $col}' >  tmpFile
	
		else
			# Find the name of the MCMC file
			name=$(head -1 ${file}|awk -F "," -v col=$i '{printf $col}')
			outfile=${name}/output

			# Finds the node number for a latent tip
			nodeNum=$(grep --text -w ${seq} ${outfile} | head -3 |tail -1 | awk '{print $2}')
	
			# Finds the column number in the mcmc that matches the "seq" latent sequence
			colNum=$(head -1 ${name}/mcmc.txt | tr '\t' '\n' | cat -n | grep t_n${nodeNum}$ |awk '{print $1}')
	
			# Put the output of the mcmc into a temporary file
			awk -v c1=$colNum '{printf "%s\n", $c1}' ${name}/mcmc.txt > tmpFile


			# Check if the column to be added is longer than columns already 
			# added to the csv file
			if [ "$i" -gt  2 ]
			then
				len1=$(wc ${seq1}.txt | awk '{print $1}')
				len2=$(wc tmpFile | awk '{print $1}')

				# If the column to be added is longer, add the appropriate number
				# of commas to the end of the csv file so the columns line up
				if [ "$len2" -gt  "$len1" ] 
				then

					numCol=$(($i-2))

					for ((j=len1;j<len2;j++));
					do
						echo $(for k in $(seq 1 ${numCol}); do printf ","; done) >> ${seq1}.txt
					done
				fi
			fi
		fi 

		paste -d "," ${seq1}.txt tmpFile > out 

		# If this is the first column, remove the "," from the file
		if [ $i == 1 ]
		then
			sed -i 's/^,//g' out 
		fi

		mv out ${seq1}.txt
	
	done

done
rm tmpFile
