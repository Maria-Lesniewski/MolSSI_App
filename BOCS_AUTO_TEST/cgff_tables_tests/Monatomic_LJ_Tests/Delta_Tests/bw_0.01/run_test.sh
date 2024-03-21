echo "Running Monotomic Delta_0.01 Test"

cgff > output.log 

echo "cgff complete, tables tool activated"

tables f_forces.Pair_Interaction.total.Pair_LJ_LJ.dat table_LJ_LJ.xvg nb delta 2.0

if [ ! -s table_LJ_LJ.xvg ]; then
        {
        echo "Failed on Table Generation"
        exit 1
        }
fi

echo "tables nonempty, comparing to standard"

sh ../../../../testing_suite_tools/compare_output2.sh ../../Data_to_Analyze/LJ_standard.xvg table_LJ_LJ.xvg 0.37 2.0 2.45 


exit_code_comparison=$?

if [[ $exit_code_comparison -eq 1 ]]; then
        {
        echo "This test failed"
        exit 1
        }
else
        echo "This test passed"
        exit 0
fi
