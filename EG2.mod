param n := 1000;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj: 0.5*sin(x[n]^2) + sum{i in 1..n-1} sin(x[1]+x[i]^2-1);
