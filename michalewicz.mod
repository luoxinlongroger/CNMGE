param pi := 3.1415926;
param n := 2;
param m := 10;
var x{i in 1..n} >=0 ,<=pi;
for{i in 1..n} 
{let x[i] := 1};

minimize obj:  -sum{i in 1..n} (sin(x[i])*sin(i*x[i]^2/pi)^(2*m));
