param n := 4;
param m := 3;
param A{1..n,1..m};
let A[1,1] := 3;
let A[2,1] := 0.1;
let A[3,1] := 3;
let A[4,1] := 0.1;
let A[1,2] := 10;
let A[2,2] := 10;
let A[3,2] := 10;
let A[4,2] := 10;
let A[1,3] := 30;
let A[2,3] := 35;
let A[3,3] := 30;
let A[4,3] := 35;

param P{1..n,1..m};
let P[1,1] := 0.3689;
let P[2,1] := 0.4699;
let P[3,1] := 0.1091;
let P[4,1] := 0.0381;
let P[1,2] := 0.117;
let P[2,2] := 0.4387;
let P[3,2] := 0.8732;
let P[4,2] := 0.5743;
let P[1,3] := 0.2673;
let P[2,3] := 0.747;
let P[3,3] := 0.5547;
let P[4,3] := 0.8828;

param alpha{i in 1..n};
let alpha[1] := 1;
let alpha[2] := 1.2;
let alpha[3] := 3.0;
let alpha[4] := 3.2;

var x{i in 1..m} >=0, <=1;
for{i in 1..m} 
{let x[i] := 1};



minimize obj: -sum{i in 1..n} (alpha[i]*exp(-sum{j in 1..m} (A[i,j]*(x[j]-P[i,j])^2)));
#
