param pi := 3.1415926;
var x1 := 2, >=-10, <=10;
var x2 := 2, >=-10, <=10;



minimize obj: -abs(sin(x1)*cos(x2)*exp(abs(1 - sqrt(x1^2+x2^2)/pi)));
#
