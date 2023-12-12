param n := 2;
var x{i in 1..n} ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj: (x[1]+2*x[2]-7)^2+(2*x[1]+x[2]-5)^2;
