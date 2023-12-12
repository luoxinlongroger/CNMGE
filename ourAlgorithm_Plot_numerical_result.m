% 消耗时间对比图
NumOfEx = 1:1:34;
%large scale problem
CnmTime_CNMGE = [1.59E+01 1.60E+02 3.67E+02 6.83E+01 4.32E+00 5.78E+01 2.84E+01 1.08E+01 2.23E-01 1.34E+01 6.00E+00 4.96E+01 5.17E+01 2.03E+02 5.52E+00 1.87E+02 7.35E-01 2.15E+00 6.29E+01 5.89E+00 4.75E+01 1.44E+01 6.83E+00 3.98E+01 7.66E+01 3.92E+02 8.85E+01 3.58E+01 1.06E+01 2.04E+01 6.20E-01 7.20E+01 1.14E+02 9.11E+01];    
CnmTime_CNMTrM = [2.16 1.85 15.53 5.66 15.57 6.87 2.03 0.58 0.73 1.02 2.33 4.45 24.83 35.06 20.11 3.22 1.87 0.76 15.10 11.53 20.02 1.62 8.05 3.14 4.69 11.03 7.78 3.92 2.03 20.52 0.89 13.34 9.06 35.75]; 
CnmTime_CNMGE_ag = [2.49E+01 7.17E+01 1.50E+02 3.03E+02 5.70E+00 8.50E+01 1.11E+02 1.16E+01 1.82E-01 2.59E+02 5.11E+00 4.76E+01 4.61E+01 1.40E+02 1.13E+01 4.12E+02 3.95E-01 2.64E+00 1.19E+02 6.02E+00 1.57E+02 1.53E+02 5.59E+00 3.38E+01 5.16E+01 2.02E+01 5.97E+01 1.63E+01 1.40E+01 3.05E+01 5.18E-01 4.94E+01 3.47E+01 2.43E+02];
CnmTime_CNMDTM = [1.21E+01 1.59E+02 3.62E+02 6.23E+02 4.72E+00 5.29E+02 2.15E+01 1.17E+01 1.13E-01 1.40E+01 6.58E+00 4.92E+01 5.33E+01 1.73E+02 5.61E+00 1.86E+02 6.05E-01 2.50E+00 6.15E+01 2.13E+00 3.14E+01 1.39E+01 3.93E+00 4.07E+01 7.37E+01 3.59E+02 8.29E+01 3.06E+01 8.76E+01 1.61E+01 3.10E-01 6.75E+01 7.42E+01 9.11E+01];
% CnmTime_MCS = [2.01E+02 1.28E+04 1.10E+02 2.60E+02 1.25E+02 3.17E+02 5.45E+01 1.26E+04 1.26E+04 1.43E+04 1.44E+02 8.48E+02 2.18E+02 2.05E+02 1.27E+04 3.34E+02 1.22E+04 1.31E+04 1.59E+02 7.88E+01 6.79E+01 7.79E+01 9.22E+01 9.52E+01 7.05E+01 1.27E+04 1.41E+02 2.17E+02 7.52E+01 2.79E+02 1.25E+04 8.24E+01 1.01E+02 1.70E+02];

figure(1)
plot(NumOfEx,CnmTime_CNMGE,'x-r');
%semilogy(NumOfEx,CnmTime_Ptc,'+-r');
hold on;
plot(NumOfEx,CnmTime_CNMTrM,'d--b');
%semilogy(NumOfEx,CnmTime_Sqp,'d-b');
hold on;
plot(NumOfEx,CnmTime_CNMDTM,'^-k');
%semilogy(NumOfEx,CnmTime_BDF,'^-k');
hold on;
plot(NumOfEx,CnmTime_CNMGE_ag,'*-g');
% semilogy(NumOfEx,CnmTime_PM,'.-g');
hold on;
% plot(NumOfEx,CnmTime_MCS,'.-c');
% % semilogy(NumOfEx,CnmTime_PM,'.-c');
% hold off
%title('The error propagations of three algorithms','Color','k','FontSize',10);
% set(gca,'XTick',2:5:41);
% set(gca,'XTickLabel',{'2','3','5','6','7','8','9','11','12','13','15','16','17','18','19','20','21','22','23','24','25','27','28','30','32','33','34','35','36','37','38','39','40','41','42','44','45','46','47'});
set(gca,'YScale','log')     

xlabel('Problem','FontSize',10);
ylabel('Computational Time','FontSize',10);
legend('CNMGE','CNMTrM','CNMDTM','CNMGE-AG')

%small scale problem

NumOfEx = 35:1:68;
CnmTime_CNMGE = [1.98E+01 3.87E-01 1.40E+00 2.38E-01 1.97E+00 3.76E-02 4.42E-02 2.66E-02 3.58E-02 3.62E-02 5.67E-01 1.86E-01 1.00E-01 2.48E+00 3.93E-01 9.51E-02 9.84E-02 3.88E+00 3.82E-01 2.88E+00 8.16E-02 5.82E-02 1.83E-01 2.15E+00 3.14E+00 1.21E+00 6.77E-01 2.80E-01 5.12E-01 7.09E-02 2.13E+00 6.80E-02 5.32E-01 6.69E-02];      
CnmTime_CNMTrM = [0.2243 0.0162 0.0255 0.0181 0.0117 0.0180 0.0140 0.0106 0.0091 0.0101 0.0141 0.0193 0.0399 0.0350 0.0185 0.0240 0.0306 0.0350 0.0493 0.0628 0.0257 0.0165 0.0438 0.0182 0.0726 0.0627 0.0351 0.0238 0.0294 0.0197 0.0230 0.0803 0.0143 0.0319 ]; 
CnmTime_CNMDTM = [1.64E+00 3.11E-01 1.24E+00 1.04E-01 1.93E+00 3.79E-02 3.45E-02 2.20E-02 1.84E-02 1.22E-02 5.59E-01 1.29E-01 8.98E-02 2.66E+00 3.80E-01 8.68E-02 8.12E-02 4.16E+00 3.68E-01 2.98E+00 6.40E-02 4.36E-02 1.19E-01 2.19E+00 3.20E+00 1.20E+00 6.77E-01 4.19E-01 4.83E-01 6.17E-02 2.25E+00 5.30E-02 5.20E-01 8.22E-02];
CnmTime_CNMGE_ag = [7.95E-02 2.18E-01 1.39E-01 7.00E-02 8.44E-01 2.04E-02 1.17E-02 1.08E-02 9.17E-03 4.41E-03 3.30E-01 2.15E-01 1.73E-02 6.04E-01 1.16E-01 1.51E-02 5.39E-02 8.18E-02 1.94E-01 1.92E+00 4.28E-02 2.17E-02 4.45E-01 2.84E+00 7.76E-02 2.13E-01 2.46E-01 6.31E-02 1.15E-01 4.67E-01 9.39E-01 3.82E-02 4.78E-01 2.66E-02];
% CnmTime_MCS = [6.38E-01 7.01E-01 3.26E-01 2.82E-01 6.13E-01 1.42E-01 4.91E-01 4.46E-01 2.48E-01 9.18E-02 2.30E-01 1.07 5.95E-01 6.62E-01 1.66E-01 1.15E-01 3.75E-01 5.65E-01 9.25E-01 3.93E-01 7.37E-01 2.57E-01 2.44E-01 7.22E-01 6.22E-01 5.86E-01 4.56E-01 4.70E-01 3.11E-01 1.59E-01 3.40E-01 4 2.64E-01 3.39E-01];

figure(2)
plot(NumOfEx,CnmTime_CNMGE,'x-r');
%semilogy(NumOfEx,CnmTime_Ptc,'+-r');
hold on;
plot(NumOfEx,CnmTime_CNMTrM,'d--b');
%semilogy(NumOfEx,CnmTime_Sqp,'d-b');
hold on;
plot(NumOfEx,CnmTime_CNMDTM,'^-k');
%semilogy(NumOfEx,CnmTime_BDF,'^-k');
hold on;
plot(NumOfEx,CnmTime_CNMGE_ag,'*-g');
% semilogy(NumOfEx,CnmTime_PM,'.-g');
hold on;
% plot(NumOfEx,CnmTime_MCS,'.-c');
% % semilogy(NumOfEx,CnmTime_PM,'.-c');
% hold off
%title('The error propagations of three algorithms','Color','k','FontSize',10);
% set(gca,'XTick',2:5:41);
% set(gca,'XTickLabel',{'2','3','5','6','7','8','9','11','12','13','15','16','17','18','19','20','21','22','23','24','25','27','28','30','32','33','34','35','36','37','38','39','40','41','42','44','45','46','47'});
set(gca,'YScale','log')     

xlabel('Problem','FontSize',10);
ylabel('Computational Time','FontSize',10);
legend('CNMGE','CNMTrM','CNMDTM','CNMGE-AG')


