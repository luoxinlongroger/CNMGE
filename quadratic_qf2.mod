param n := 1000;
var x{i in 1..n} >=0;
for{i in 1..n} 
{let x[i] := 2};



minimize obj: - x[n] + 0.5* sum{i in 1..n}(i*(x[i]^2-1)^2);
