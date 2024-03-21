# Author Maria Lesniewski
# This is a script designed to determine if 2 tabulated potentials are the same within a passed tolerance within a particular region
# Usage : bash compare_output.sh <File_1> <File_2> <x_domainb> <x_domaine> <tolerance>
# Explain: File_X are the files to compare to eachother, x_domain b and e are the begining and end of the domains for which you expect
#          The files compared to reasonably agree, tolerance is the pointwise error that if achieved will cause this script to exit 
# This script is specifically designed to compare the output of the tables tool inside of BOCS when analytic and tabulate potentials are to be defined along the same bin spacing
# This script could also probably be made more efficient by piping, but I have tried to maximize readability 
#!/bin/bash

#Check if args passed:
if [ "$#" -lt 5 ] || [ ${1} == "-h" ]
then
        echo "Too few arguments passed or -h, did not execute"
        echo "Pass 1.File1 2.File2 3. x_1 4. x_2 5. tolerance"
	echo "Files to be compared, domain, pointwise error"
        exit 0
fi

file1=${1}
file2=${2}
begin=${3}
end=${4}
tol=${5}

# First We need to trim the domains of File1 and File2 to be defined along the same 
awk -v a="$begin" -v b="$end" '$1 >= a && $1 <= b {print $1 " " $6 " " $7}' "$file1" > file_1_trimmed.txt
awk -v a="$begin" -v b="$end" '$1 >= a && $1 <= b {print $1 " " $6 " " $7}' "$file2" > file_2_trimmed.txt

# Then We need to compare line by line the difference in forces and potentials
#Checking the potentials first
awk 'FNR==NR { y1[FNR] = $2; next} { y2[FNR] = $2} END { for (i = 1; i <=FNR; i++) print y1[i] - y2[i] }' file_1_trimmed.txt file_2_trimmed.txt > dif_potential.txt
awk 'FNR==NR { y1[FNR] = $3; next} { y2[FNR] = $3} END { for (i = 1; i <=FNR; i++) print y1[i] - y2[i] }' file_1_trimmed.txt file_2_trimmed.txt > dif_forces.txt

# This is a pointwise error check, so we'll find the maximum error and determine if it is out of spec 
max_potential_error=$(awk 'NR == 1 || sqrt($1^2) > max{max=sqrt($1^2)} END {print max}' dif_potential.txt)
max_force_error=$(awk 'NR == 1 || sqrt($1^2) > max{max=sqrt($1^2)} END {print max}' dif_forces.txt)

if [[ $max_force_error > $tol ]]; then
	{
	echo "ERROR POINTWISE FORCES WERE NOT WITHIN TOLERANCE: $tol"
	echo "Max FORCE ERROR: $max_force_error"
	echo "Max POTENTIAL ERROR: $max_potential_error"
	exit 1
	}	
else
	echo "PASSED"
	exit 0
fi

