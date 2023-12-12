param n := 1000;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj: sum{i in 1..n/2} (((x[2*i-1]-3)/100)^2 - (x[2*i-1] - x[2*i]) + exp(20*(x[2*i-1]-x[2*i])));
