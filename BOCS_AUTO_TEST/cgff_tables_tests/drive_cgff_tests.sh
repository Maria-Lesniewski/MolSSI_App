#!/bin/bash

echo "LAUNCHING CGFF TESTS DEFAULT TESTING SYTEMS 1. (a-f) 2. 3. "
# LAUNCH 1 #######################################################################################
# LAUNCH 1a ######################################################################################
cd Monatomic_LJ_Tests/
echo "1. LJ Monatomic Systems [NPT]"
cd Bspline_Test/
echo "a. Bspline Test"
sh run_test.sh > Bspline_testoutput.txt
test_exit=$?
if [[ $test_exit -eq 1 ]]; then
        {
        echo "FAILED"
	bresult=1
        }
else
        echo "PASSED"
	bresult=0
fi

#LAUNCH 1b #####################################################################################
cd ../Delta_Tests/bw_0.01/
echo "b. Delta 0.01 Test" 
sh run_test.sh > Delta_testoutput.txt
test_exit=$?
if [[ $test_exit -eq 1 ]]; then
        {
        echo "FAILED"
        D01=1
        }
else
        echo "PASSED"
        D01=0
fi
#LAUNCH 1c #####################################################################################
echo "c. Delta 0.004 Test"
cd ../bw_0.004/
sh run_test.sh > Delta_testoutput.txt
test_exit=$?
if [[ $test_exit -eq 1 ]]; then
        {
        echo "FAILED"
        D001=1
        }
else
        echo "PASSED"
        D001=0
fi
#LAUNCH 1d #####################################################################################
cd ../../Linear_Test/bw_0.01/
echo "d. Linear 0.01 Test"
sh run_test.sh > linear_testoutput.txt
test_exit=$?
if [[ $test_exit -eq 1 ]]; then
        {
        echo "FAILED"
        L01=1
        }
else
        echo "PASSED"
        L01=0
fi
#LAUNCH 1e #####################################################################################
echo "e. Linear 0.003 Test"
cd ../bw_0.003/
sh run_test.sh > linear_testoutput.txt
test_exit=$?
if [[ $test_exit -eq 1 ]]; then
        {
        echo "FAILED"
        L001=1
        }
else
        echo "PASSED"
        L001=0
fi
#LAUNCH 1f #####################################################################################
echo "f. Power Test"
cd ../../Power_Test/
sh run_test.sh > Power_testoutput.txt
test_exit=$?
if [[ $test_exit -eq 1 ]]; then
        {
        echo "FAILED"
        power=1
        }
else
        echo "PASSED"
        power=0
fi
###########################################SUMMARY############################################
echo "cgff_tables_test_summary"
echo "0 = passed, 1 = Failed"
echo "Bspline     : $bresult"
echo "Delta 0.01  : $D01"
echo "Delta 0.001 : $D001"
echo "Linear 0.01 : $L01"
echo "Linear 0.001: $L001"
echo "Power       : $power"
