#!/bin/bash

sed -i "s:linear:Bspline:g" par.txt
sed -i "s:0.500    0.000      180        0:3.000    0.000      180        4   0:g" par.txt #angles
sed -i "s:0.500     -180      180        0:3.000     -180      180        4   0:g" par.txt #Dihedrals
sed -i "s:0.002    0.000      0.7        0:0.02    0.000      0.7         4   0:g" par.txt #Bonds
