param pi := 3.1415926;
var x1 := 2, >=-100, <=100;
var x2 := 2, >=-100, <=100;



minimize obj:  x1^2+2*x2^2-0.3*cos(3*pi*x1)-0.4*cos(4*pi*x2)+0.7 ;
#
