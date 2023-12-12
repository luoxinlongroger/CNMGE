param n := 10;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj:   sum{j in 1..n} x[j]^2;