clear;
close;

%% PART I- DATA IMPORT

%% Descret DATA

% %% Initialize variables.
% 
% % For lab computer
% 
% path = 'C:\Documents and Settings\nilesh\My Documents\MATLAB\MKY1';
%                 % filename = 'path''C:\Documents and Settings\nilesh\My Documents\MATLAB\MKY1\modelfit.txt';
% filename = strcat(path,'\SH Programs\TIV\TIV_Estrin\equispaced.txt');
% 
% % for home computer
% % filename = 'C:\Users\Intel\Documents\MATLAB\TIV Feb 19\realTIV\equispaced.txt'; %equispaced modelfit.txt
% 
% %%
% 
% delimiter = '\t';
% formatSpec = '%f%f%[^\n\r]'; 
% 
% %% Open the text file.
% fileID = fopen(filename,'r');
% 
% %% Read columns of data according to format string.
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
% 
% %% Close the text file.
% fclose(fileID);
% 
% 
% %% Create output variable
% 
% data1 = dataset(dataArray{1:end-1}, 'VarNames', {'trueStrain','trueStress'});
% data = iddata(data1.trueStress(1:100), [], 0.0022,'Name', 'true stress'); % Ts = 0.0022 
% 
% %data = iddata(data1.trueStress(1:100),data1.trueStrain(1:100)); % [o/p, i/p]
% 
% %% Clear temporary variables
% 
% clearvars filename delimiter formatSpec fileID dataArray ans;
% set(data, 'OutputName', 'Stress');
% set(data, 'OutputUnit', 'MPa');
% set(data, 'Tstart',0.0024, 'TimeUnit', 's');%0.0024

%% Continulous data

% For lab computer

path = 'C:\Documents and Settings\nilesh\My Documents\MATLAB\MKY1';
                % filename = 'path''C:\Documents and Settings\nilesh\My Documents\MATLAB\MKY1\modelfit.txt';
filename = strcat(path,'\SH Programs\solutionized.txt'); %2hr aged

% for home computer
% filename = 'C:\Users\Intel\Documents\MATLAB\TIV Feb 19\realTIV\equispaced.txt'; %equispaced modelfit.txt

g = fopen(filename,'r');                            

dt = textscan(g, '%f %f');

strain = dt{1,1}(:);
stress = dt{1,2}(:);
w = polyfit(strain,stress,6);
t =  0.0:0.002:0.224; %0.0:0.002:0.226; %0.0024:0.0002:0.2264;%
Stress = polyval(w,t');
data = iddata(Stress,[],0.002,'Name','truestress');


%% PART - II MODEL


%% Parameter definition

parameters = {1,1,1,1,1,1,1,1,1,1,1};

%% Order of the model % Order 
% Vector with three entries [Ny Nu Nx], specifying the number of model
% outputs Ny, the number of inputs Nu, and the number of states Nx

order         = [1 0 2];  
                           
                           
%% State variable initialization


rho_m0 = 1e11;
rho_f0 = 1e10; 

initial_states = [rho_m0; rho_f0];

%% Model definition


TIVmodel = idnlgrey('TIVSum', order, parameters, initial_states);   % for continuous data

TIVmodel.Algorithm.MaxIter = 100;


%%  Fixiing of parameters             

setpar(TIVmodel, 'Fixed', {false false false false false true true true true true true});

%% Fixing of initial state variables

setinit(TIVmodel,'Fixed',{true true});

%% To fix algorithm options. If we dont fix then it uses auto mode to select the algorithm options

 set(TIVmodel,'TimeUnit','s');
 set(TIVmodel,'CovarianceMatrix','Estimate');
 TIVmodel.Algorithm.SearchMethod = 'lm'; %'lm',gna,gn,grad,auto
 
%  lavenberg-Marquardt method setting

% lm.StartValue = 0.001;
% lm.LMStep = 10;
% lm.MaxBisections = 25;
% lm.RelImprovement = 1e-9;

               
TIVmodel.Algorithm.SimulationOptions.Solver = 'ode45';
%'ode45' 'ode23' 'ode113' 'ode15s' 'ode23s' 'ode23tb' 'ode5' 'ode4' 'ode3' ode2' 'ode1' 'auto'

% TIVmodel.Algorithm.GradientOptions.DiffScheme = 'Central approximation';

%% Setting parameters' values and limits


TIVmodel.Parameters(1).Value = 7e7;%5e7;             %k1
TIVmodel.Parameters(2).Value = 8;%1.73;              %k2
TIVmodel.Parameters(3).Value = 600;%11.4;              %k3
TIVmodel.Parameters(4).Value = 46;%53;             %kL
TIVmodel.Parameters(5).Value= 6e-6;             %Ls
TIVmodel.Parameters(6).Value = 30e-6;           %L0
TIVmodel.Parameters(7).Value = 2.96;            %M
TIVmodel.Parameters(8).Value = 2.86e-10;        %b
TIVmodel.Parameters(9).Value = 0.33;%0.2;             %alpha
TIVmodel.Parameters(10).Value = 89; % 93.86; 262         %sigma_i
TIVmodel.Parameters(11).Value = 27e3;           %G MPa

%  estimated value of k1 is 70000000.001053 +- 74875113.213859 
%  estimated value of k2 is 8.103777 +- 5.628366 
%  estimated value of k3 is 559.041588 +- 780.334093 
%  estimated value of kL is 46.722046 +- 27.511561 
%  estimated value of Ls is 0.000018 +- 0.000005 


% TIVmodel.Parameters(1).Minimum = 1e5;
% TIVmodel.Parameters(2).Minimum = 1;
% TIVmodel.Parameters(3).Minimum = 1;
% TIVmodel.Parameters(4).Minimum = 1;
% TIVmodel.Parameters(5).Minimum = 0;

% TIVmodel.Parameters(1).Maximum = 0;
% TIVmodel.Parameters(2).Maximum = 0;
% TIVmodel.Parameters(3).Maximum = 0;
% TIVmodel.Parameters(4).Maximum = 0;
% TIVmodel.Parameters(5).Maximum = 0;

%% Setting algorithm properties

%TIVModel.Algorithm.Advance.GnPinvConst = 1e2;



%% Estimation of the model using pem function

TIVmodel = pem(data,TIVmodel,'Display','Full'); 


%% Outputs of the model

% b_est = TIVmodel.Parameters(4).Value;
[nonlinear_b_est, sdnonlinear_b_est] = ...
                            getpvec(TIVmodel, 'free');

B = fopen('SHmodel_result.txt','a');
fprintf(B,'Solver: %s \n\n',TIVmodel.Algorithm.SimulationOptions.Solver);
fprintf(B,' estimated value of k1 is %f +- %f \n', nonlinear_b_est(1), sdnonlinear_b_est(1));
fprintf(B,' estimated value of k2 is %f +- %f \n', nonlinear_b_est(2), sdnonlinear_b_est(2));
fprintf(B,' estimated value of k3 is %f +- %f \n', nonlinear_b_est(3), sdnonlinear_b_est(3));
fprintf(B,' estimated value of kL is %f +- %f \n', nonlinear_b_est(4), sdnonlinear_b_est(4));
fprintf(B,' estimated value of Ls is %f +- %f \n\n', nonlinear_b_est(5), sdnonlinear_b_est(5));

fclose(B);

                        
%% Plotting of the results

compare(data,TIVmodel,'-*b')

[y1,fit,x0]= compare(data,TIVmodel); 
fit1 = goodnessOfFit(Stress,y1.y,'NRMSE'); %MSE,NMSE

%% Some commands which can be given at the prompt to extract more results

%  sim(TIVmodel,data)
% findstates(TIVmodel,data,initial_states)
% sys = n4sid(data,4) 

%% Solver options

% ODE (Ordinary Differential/Difference Equation) solver
% for solving state space equations.A. Variable-step

% solvers for time-continuous idnlgrey models:'ode45' — Runge-Kutta (4,5)
% solver for nonstiff problems.'ode23' — Runge-Kutta (2,3)
% solver for nonstiff problems.'ode113' — Adams-Bashforth-Moulton
% solver for nonstiff problems.'ode15s' — Numerical Differential
% Formula solver for stiff problems.'ode23s' — Modified Rosenbrock
% solver for stiff problems.'ode23t' — Trapezoidal solver
% for moderately stiff problems.'ode23tb' — Implicit Runge-Kutta
% solver for stiff problems.B. Fixed-step solvers for time-continuous idnlgrey models:'ode5' — Dormand-Prince
% solver.'ode4' — Fourth-order Runge-Kutta
% solver.'ode3' — Bogacki-Shampine
% solver.'ode2' — Heun or improved
% Euler solver.'ode1' — Euler solver.C. Fixed-step solvers for time-discrete idnlgrey models: 'FixedStepDiscrete'D.
% General: 'Auto' — Automatically chooses
% one of the previous solvers (default).
                        
