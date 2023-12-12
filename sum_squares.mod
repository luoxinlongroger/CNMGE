param n := 1000;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj:   sum{j in 1..n} j*x[j]^2;