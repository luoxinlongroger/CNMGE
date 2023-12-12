function [obj, grad] = generatedObjective(inputVariables)
%generatedObjective Compute objective function value and gradient
%
%   OBJ = generatedObjective(INPUTVARIABLES) computes the objective value
%   OBJ at the point INPUTVARIABLES.
%
%   [OBJ, GRAD] = generatedObjective(INPUTVARIABLES) additionally computes
%   the objective gradient value GRAD at the current point.
%
%   Auto-generated by prob2struct on 06-Oct-2022 14:42:18

%% Map solver-based variables to problem-based.
x = inputVariables(:);

%% Compute objective function.
arg2 = 14;
arg3 = x(1);
arg5 = x(1);
arg6 = 3;
arg7 = arg5.^2;
arg9 = 14;
arg10 = x(2);
arg12 = 6;
arg13 = x(1);
arg15 = (arg12 .* arg13);
arg16 = x(2);
arg18 = x(2);
arg19 = 3;
arg20 = arg18.^2;
arg23 = ((x(1) + x(2)) + 1);
arg24 = arg23.^2;
arg25 = (((((19 - (arg2 .* arg3)) + (arg6 .* arg7)) - (arg9 .* arg10)) + (arg15 .* arg16)) + (arg19 .* arg20));
arg27 = 32;
arg28 = x(1);
arg30 = x(1);
arg31 = 12;
arg32 = arg30.^2;
arg34 = 48;
arg35 = x(2);
arg37 = 36;
arg38 = x(1);
arg40 = (arg37 .* arg38);
arg41 = x(2);
arg43 = x(2);
arg44 = 27;
arg45 = arg43.^2;
arg47 = 2;
arg48 = x(1);
arg50 = 3;
arg51 = x(2);
arg52 = ((arg47 .* arg48) - (arg50 .* arg51));
arg53 = arg52.^2;
arg54 = (((((18 - (arg27 .* arg28)) + (arg31 .* arg32)) + (arg34 .* arg35)) - (arg40 .* arg41)) + (arg44 .* arg45));
arg55 = (1 + (arg24 .* arg25));
arg56 = (30 + (arg53 .* arg54));
obj = (arg55 .* arg56);

if nargout > 1
    %% Compute objective gradient.
    % To call the gradient code, notify the solver by setting the
    % SpecifyObjectiveGradient option to true.
    arg57 = (1.*arg55(:));
    arg58 = ((-((arg57.*arg54(:)).*2.*(arg52(:)))).*arg50(:));
    arg59 = zeros([2,size(arg58,2)]);
    arg59(2,:) = arg58;
    arg60 = (((arg57.*arg54(:)).*2.*(arg52(:))).*arg47(:));
    arg61 = zeros([2,size(arg60,2)]);
    arg61(1,:) = arg60;
    arg62 = (arg57.*arg53(:));
    arg63 = ((arg62.*arg44(:)).*2.*(arg43(:)));
    arg64 = zeros([2,size(arg63,2)]);
    arg64(2,:) = arg63;
    arg65 = ((-arg62).*arg40(:));
    arg66 = zeros([2,size(arg65,2)]);
    arg66(2,:) = arg65;
    arg67 = (((-arg62).*arg41(:)).*arg37(:));
    arg68 = zeros([2,size(arg67,2)]);
    arg68(1,:) = arg67;
    arg69 = arg62;
    arg70 = (arg69.*arg34(:));
    arg71 = zeros([2,size(arg70,2)]);
    arg71(2,:) = arg70;
    arg72 = arg69;
    arg73 = ((arg72.*arg31(:)).*2.*(arg30(:)));
    arg74 = zeros([2,size(arg73,2)]);
    arg74(1,:) = arg73;
    arg75 = ((-arg72).*arg27(:));
    arg76 = zeros([2,size(arg75,2)]);
    arg76(1,:) = arg75;
    arg77 = (1.*arg56(:));
    arg78 = ((arg77.*arg25(:)).*2.*(arg23(:)));
    arg79 = arg78;
    arg80 = arg79;
    arg81 = zeros([2,size(arg80,2)]);
    arg81(2,:) = arg80;
    arg82 = arg79;
    arg83 = zeros([2,size(arg82,2)]);
    arg83(1,:) = arg82;
    arg84 = (arg77.*arg24(:));
    arg85 = ((arg84.*arg19(:)).*2.*(arg18(:)));
    arg86 = zeros([2,size(arg85,2)]);
    arg86(2,:) = arg85;
    arg87 = arg84;
    arg88 = (arg87.*arg15(:));
    arg89 = zeros([2,size(arg88,2)]);
    arg89(2,:) = arg88;
    arg90 = ((arg87.*arg16(:)).*arg12(:));
    arg91 = zeros([2,size(arg90,2)]);
    arg91(1,:) = arg90;
    arg92 = ((-arg87).*arg9(:));
    arg93 = zeros([2,size(arg92,2)]);
    arg93(2,:) = arg92;
    arg94 = arg87;
    arg95 = ((arg94.*arg6(:)).*2.*(arg5(:)));
    arg96 = zeros([2,size(arg95,2)]);
    arg96(1,:) = arg95;
    arg97 = ((-arg94).*arg2(:));
    arg98 = zeros([2,size(arg97,2)]);
    arg98(1,:) = arg97;
    grad = arg59 + arg61 + arg64 + arg66 + arg68 + arg71 + arg74 + arg76 + arg81 + arg83 + arg86 + arg89 + arg91 + arg93 + arg96 + arg98;
    
end
end