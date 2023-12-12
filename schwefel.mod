param n := 1000;
var x{i in 1..n} >=-500, <=500;
for{i in 1..n} 
{let x[i] := 1};



minimize obj: 418.9829*n + sum{j in 1..n} (x[j]*sin(sqrt(sqrt(x[j]^2))));
# 