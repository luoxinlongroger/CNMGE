param n := 1000;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj: sum{i in 1..n/2} ((x[2*i-1]^2+x[2*i]^2+x[2*i-1]*x[2*i])^2 + sin(x[2*i-1])^2 + cos(x[2*i])^2);
