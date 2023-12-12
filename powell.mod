param n := 1000;
var x{i in 1..n} >=-4, <=5;
for{i in 1..n} 
{let x[i] := 2};


minimize obj:   sum{i in 1..n/4} ((x[4*i-3]+10*x[4*i-2])^2 + 5*(x[4*i-1]-x[4*i])^2 + (x[4*i-2]-2*x[4*i-1])^4 + 10*(x[4*i-3]-x[4*i])^4);