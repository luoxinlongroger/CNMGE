param n := 1000;
var x{i in 1..n} >=-5,<=5;
for{i in 1..n} 
{let x[i] := 2};



minimize obj:   0.5*(sum{j in 1..n} x[j]^4) - 0.5*(16*sum{j in 1..n} x[j]^2)+ 0.5*(5*sum{j in 1..n}x[j]);