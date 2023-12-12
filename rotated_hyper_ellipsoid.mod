param n := 1000;
var x{i in 1..n} >=-65.536,<=65.536;
for{i in 1..n} 
{let x[i] := 2};



minimize obj:    sum{j in 1..n}(sum{i in 1..j} x[i]^2);