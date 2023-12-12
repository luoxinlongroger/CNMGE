param n := 1000;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj:   (x[1]-1)^2 + sum{j in 2..n} j*(2*x[j]^2-x[j-1])^2 ;