param n := 4;
param beta := 4;
var x{i in 1..n} >=-4, <=4;
for{i in 1..n} 
{let x[i] := 1};



minimize obj:  sum{i in 1..n} ((sum{j in 1..n} (j^i+beta)*((x[j]/j)^i-1))^2);
#
