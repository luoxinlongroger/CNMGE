param pi := 3.1415926;
param t := 1 / (8*pi);
param s := 10;
param r := 6;
param c := 5/pi;
param b := 5.1 / (4*pi^2);
param a := 1;

var x1 := 2;
var x2 := 2;


minimize obj:  a * (x2 - b*x1^2 + c*x1 - r)^2 + s*(1-t)*cos(x1) + s;
