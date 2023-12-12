param n := 1000;
param noise := 0.5;
var x{i in 1..n};
for{i in 1..n} 
{let x[i] := 2};


minimize obj:   noise + sum{i in 1..n} x[i]^4;