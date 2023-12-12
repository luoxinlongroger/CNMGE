
var x1 := 2;
var x2 := 2;
var x3 := 2;


minimize obj:  sum{i in 1..3} (exp(-0.1*i*x1)-exp(-0.1*i*x2)-(exp(-0.1*i)-exp(-i))*x3)^2;