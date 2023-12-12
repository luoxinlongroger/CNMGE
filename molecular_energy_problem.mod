param n := 1000;
param k := 4.141720682;
param r := 10.60099896;
var x{i in 1..n} >=0;
for{i in 1..n} 
{let x[i] := 1};



minimize obj:  1000 + sum{j in 1..n} ((-1)^j)/sqrt(r-k*cos(x[j]))+ sum{j in 1..n} cos(3*x[j]);
# 