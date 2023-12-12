param n := 1000;
var x{i in 1..n} >=-n^2,<=n^2;
for{i in 1..n} 
{let x[i] := 2};



minimize obj:    x[1]^2 + 2*x[1] + 1000 + (sum{j in 2..n} x[j]^2)-(sum{j in 2..n} x[j]*x[j-1]) - (sum{j in 2..n} 2*x[j]);