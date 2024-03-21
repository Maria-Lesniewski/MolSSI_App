# Author Maria Lesniewski
# This is a script designed to determine if 2 tabulated potentials are the same within a passed tolerance within a particular region
# Usage : bash compare_output.sh <File_1> <File_2> <x_domainb> <x_domaine> <tolerance>
# Explain: File_X are the files to compare to eachother, x_domain b and e are the begining and end of the domains for which you expect
#          The files compared to reasonably agree, tolerance is the pointwise error that if achieved will cause this script to exit 
# This script is specifically designed to compare the output of the tables tool inside of BOCS when analytic and tabulate potentials are passed in a partic order
# We assume that the analytic standard has a much denser grid that encompasses the BOCS output 
# This script could also probably be made more efficient by piping instead of writing and removing scripts between steps, but I have tried to maximize readability 
#!/bin/bash

#Check if args passed:
if [ "$#" -lt 5 ] || [ ${1} == "-h" ]
then
        echo "Too few arguments passed or -h, did not execute"
	echo "Pass 1.File1(standard) 2.File2(test) 3. x_1 4. x_2 5. tolerance"
	echo "Files to be compared, domain, pointwise error"
        exit 0
fi

file1=${1}
file2=${2}
begin=${3}
end=${4}
tol=${5}

# First We need to trim the domains of File1 and File2 to be defined along the same well sampled region
awk -v a="$begin" -v b="$end" '$1 >= a && $1 <= b {print $1 " " $6 " " $7}' "$file1" > standard_trimmed.txt
awk -v a="$begin" -v b="$end" '$1 >= a && $1 <= b {print $1 " " $6 " " $7}' "$file2" > test_trimmed.txt

#After we've defined the domains that are well sampled, we need to account possible differences in bin spacing between the standard and test data_sets 
awk 'FNR == NR {smalldata[$1] = $2 ; next} { for (i in smalldata) if (i == $1) print $1 " " $2 " " $3 }' test_trimmed.txt standard_trimmed.txt > standard_on_test_domain.txt

# Then We need to compare line by line the difference in forces and potentials

	#Checking the potentials first
	awk 'FNR==NR { y1[FNR] = $2; next} { y2[FNR] = $2} END { for (i = 1; i <=FNR; i++) print y1[i] - y2[i] }' standard_on_test_domain.txt test_trimmed.txt > dif_potential.txt
	awk 'FNR==NR { y1[FNR] = $3; next} { y2[FNR] = $3} END { for (i = 1; i <=FNR; i++) print y1[i] - y2[i] }' standard_on_test_domain.txt test_trimmed.txt > dif_forces.txt
	#awk 'FNR==NR { standard_f[NR] = $3; next}  {print $3 - standard_f[FNR]}' standard_on_test_domain.txt test_trimmed.txt > dif_forces.txt2 # This was me sanity checking the for notation
	
	# This is a pointwise error check, so we'll find the maximum error and determine if it is out of spec 
	max_potential_error=$(awk 'NR == 1 || sqrt($1^2) > max{max=sqrt($1^2)} END {print max}' dif_potential.txt)
	max_force_error=$(awk 'NR == 1 || sqrt($1^2) > max{max=sqrt($1^2)} END {print max}' dif_forces.txt)

# Now we'll cleanup our mess
rm standard_on_test_domain.txt test_trimmed.txt standard_trimmed.txt dif_potential.txt dif_forces.txt

# Then we'll return the proper error code so that the driver script can do its thing
# If the error is zero, probably the tables tool didnt execute
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

