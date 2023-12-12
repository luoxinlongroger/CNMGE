param n := 1000;
param pi := 3.1415926;

var x{i in 1..n} >=-32.768, <=32.768;
for{i in 1..n} 
{let x[i] := 1};


minimize obj:  sin(pi*(1 + (x[1]-1)/4))^2 + sum{j in 1..n-1} ((((x[j]-1)/4)^2)*10*sin(pi*(1+(x[j]-1)/4)+1)^2) + (((x[n]-1)/4)^2)*(1+sin(2*pi*(1 + (x[n]-1)/4))^2);
# 