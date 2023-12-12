var x1 := 2, >=-10, <=10;
var x2 := 2, >=-10, <=10;


minimize obj:  (sum{i in 1..5} i*cos((i+1)*x2+i))*(sum{i in 1..5} i*cos((i-1)*x1+i));
#
