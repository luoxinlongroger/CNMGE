param n := 1000;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};


minimize obj:   sum{i in 1..n-1} (100*(x[i+1]-x[i]^2)^2 + (x[i]+1)^2);