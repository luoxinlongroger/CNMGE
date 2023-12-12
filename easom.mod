param pi := 3.1415926;
var x1 := 2, x>=-100, x<=100;
var x2 := 2, x>=-100, x<=100;


minimize obj:  -cos(x1)*cos(x2)*exp(-(x1-pi)^2-(x2-pi)^2);
#
