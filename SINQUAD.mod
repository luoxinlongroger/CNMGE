param n := 1000;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj: (x[1]-1)^4 + (x[n]^2-x[1]^2)^2 +sum{i in 2..n-1} (sin(x[i]-x[n])-x[1]^2+x[i]^2)^2;