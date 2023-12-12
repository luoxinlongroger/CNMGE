param n := 1000;
param m := n/2;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};


minimize obj:   sum{i in 1..m} ((x[2*i-1]-10)^2 + (x[2*i-1]*x[2*i]-50000)^2);