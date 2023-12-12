
Dear Colleguages,

you are welcome!

---------------------------------

The test environment has been set up. We have put our software and the test problems in the 
directory "/CNMGE_code". It has 8 subdirectories and we describes 
them as follows:

1. The subdirectory "cnmge_programs" includes the programs of our method (CNMGE) and 
it needs to run in MATLAB2020b or after this version.  The programs are put in its subdirectory 
"cnmge_programs". 

The "cnmge_programs" subdirectory includes 13 files. The main function of our method is 
CNMGE.m and it has the following structure:

"function [f_opt,x_opt,CPU_time,Num_stp] = CNMGE(func,x0,ub,lb)".

Its first input argument is the test function. The first output argument is the minimum value 
computed by CNMGE and the third output argument is its computational time.

The aim of the function CNMDTM.m (Algorithm 3 in reference [1] is to find multiple stationry 
points of the objective function from multi-start points, which is called by CNMGE module.

The aim of the function QGE.m (Algorithm 4 in reference [1]) is the implementation of the 
quasi-genetic evolution algorithm, which is called by CNMGE module.

The main of the function CNMTr (Algorithm 1 in reference [1]) is to find a stationary stationry 
point of the objective function, which is called by the CNMDTM module. 

The main of the function CNMDT (Algorithm 2 in reference [1]) is to find a new stationary stationary 
point of the objective function, which is called by the CNMDTM moule. 

These five modules are the main modules in the subdirectory "cnmge_programs". The rest two 
modules (jaco_deflation and jacobianfun) are the assistant functions and they compute the Jacobian 
matrix of the gradient, which are called by CNMDT and CNMTr, respectively. 


The CNMGE method can be tested via typing "test_CNMGE.m"  or "test_CNMGE_AG.m" in the command 
line under the MATLAB environment (MATLAB2021b) and inputting the number from 1 to 68,  which 
represents the test problem number. test_CNMGE is to test CNMGE by by  automatic differentiation 
of the reverse mode to compute the gradient of the objective function. And test_CNMGE_AG is to 
test CNMGE by the analytical gradient of the objective function, where these analytical gradients and 
the objective functions are all in the subdirectory "test_problems".
 
"test_CNMTrM.m" and "test_CNMDTM" are the other two test files and they are to test the performance 
of CNMTr with multi-start points (CNMTrM) and that of CNMDTM for global optimization problems, 
respectively. Their aims are to test the performance improvement of CNMTr with multi-start points 
(CNMTrM), CNMDTM, and (CNMDTM with GQE), respectively. 

2. The subdirectory "CMA-ES" includes the implemented files of the CMA-ES method and the test 
problems are put in "test_function". The CMA-ES method can be tested via typing "test_cmaes_main.m" 
in the command line under the MATLAB environment and inputting the number from 1 to 68, which 
represents the test problem number.

3. The subdirectory "MCS" includes the implemented files of the MCS method and the test problem 
files are put in its subdirectory "jones". The MCS method can be test via typing "test_mcs_main.m" in 
the command line under the MATLAB environment and inputting the number from 1 to 68, which 
represents the test problem number.

4. The subdirectory "Computational results of Couenne with NEOS" includes the numerical results of 
the Couenne method under the AMPL environment with the NEOS server, which is available at
https://neos-server.org/neos/. The 68 test problem files are put in the subdirectory "Couenne_mod", 
which is coded by AMPL.

6. The subdirectory "glods_0.3" includes the implemented files of the GLODS method and the test 
problems are also put in the subdirectory "test_problems". You can test GLODS just like testing 
CNMGE via typing "test_glods.m".

7. The subdirectory "VRBBO" includes the implemented files of the VRBBO method and the test 
problems are also put in the subdirectory "test_problems". You can test VRBBO just like testing 
CNMGE via typing "test_VRBBO.m".

8. The built-in subroutine GlobalSearch of the MATLAB environment can be tested via typing 
"test_globalsearch_main.m" in the command line and inputting number from 1 to 68, which represents 
the test problem number. The test problems are also put in the subdirectory "test_problems".

[1] The detailed descriptions of their design principles and the numerical results can reference to the 
manuscript entitled "Continuation Newton methods with deflation techniques and quasi-genetic 
evolution for global optimization problems, which is available at http://arxiv.org/abs/2107.13864". 
We also put its copy in the subdirectory "cnmge_code" and its file name is 
"Luo-Xiao-Zhang-arXiv-2023-Dec-11-Continuation Newton methods for global optimization.pdf". 

We hope that you have a pleasant experience for testing our software.


