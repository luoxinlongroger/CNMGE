param n := 1000;
param pi := 3.1415926;
var x{i in 1..n} >=-5.12,<=5.12;
for{i in 1..n} 
{let x[i] := 2};



minimize obj:  10*n + sum{j in 1..n} (x[j]^2-10*cos(2*pi*x[j]));