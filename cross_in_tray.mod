param pi := 3.1415926;
var x1 := 2,>=-10,<=10;
var x2 := 2,>=-10,<=10;


minimize obj:  -0.0001 * (abs(sin(x1)*sin(x2)*exp(abs(100 - sqrt(x1^2+x2^2)/pi)))+1)^0.1;