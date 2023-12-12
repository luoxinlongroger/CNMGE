param n := 10;
var x{i in 1..n} >= -600, <= 600 ;
for{i in 1..n} 
{let x[i] := 2};



minimize obj: (sum{i in 1..n} x[i]^2)/4000 -cos(x[1]/sqrt(1))*cos(x[2]/sqrt(2))*cos(x[3]/sqrt(3))*cos(x[4]/sqrt(4))*cos(x[5]/sqrt(5))*cos(x[6]/sqrt(6))*cos(x[7]/sqrt(7))*cos(x[8]/sqrt(8))*cos(x[9]/sqrt(9))*cos(x[10]/sqrt(10)) + 1;
