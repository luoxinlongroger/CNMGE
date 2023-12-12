param n := 1000;
param a := 20;
param b := 0.2;
param pi := 3.1415926;
param c := 2*pi;
var x{i in 1..n} >=-32.768, <=32.768;
for{i in 1..n} 
{let x[i] := 1};



minimize obj:  -a*exp(-b*sqrt((sum{j in 1..n}(x[j]^2))/n))-exp((sum{j in 1..n}(cos(c*x[j])))/n) + a + exp(1);
# 