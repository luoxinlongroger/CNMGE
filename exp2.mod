var x1 := 2, >=0, <=20;
var x2 := 2, >=0, <=20;



minimize obj:  sum{i in 1..10} ((exp(-(i-1)*x1/10)-5*exp(-(i-1)*(x2/10))-exp(-(i-1)/10)+5*exp(-(i-1)))^2);
#
