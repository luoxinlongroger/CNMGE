param n := 1000;
var x{i in 1..n} >=-10, <=10;
for{i in 1..n} 
{let x[i] := 2};


minimize obj:  -sum{i in 1..n} (sin(2*x[i]+1)+2*sin(3*x[i]+2)+3*sin(4*x[i]+3)+4*sin(5*x[i]+4)+5*sin(6*x[i]+5));