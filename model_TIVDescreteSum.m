clear;
close;

%% PART I- DATA IMPORT

%% Descret DATA

%% Initialize variables.

% For lab computer

path = 'C:\Documents and Settings\nilesh\My Documents\MATLAB\MKY1';
filename = strcat(path,'\SH Programs\solutionized.txt'); %2hr aged

% for home computer
% filename = 'C:\Users\Intel\Documents\MATLAB\TIV Feb 19\realTIV\equispaced.txt'; %equispaced modelfit.txt

%%

delimiter = '\t';
formatSpec = '%f%f%[^\n\r]'; 

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);


%% Create output variable

data1 = dataset(dataArray{1:end-1}, 'VarNames', {'trueStrain','trueStress'});
%data = iddata(data1.trueStress(1:100), [], 0.0022,'Name', 'true stress'); % Ts = 0.0022 

data = iddata(data1.trueStress(1:100),data1.trueStrain(1:100)); % [o/p, i/p]

%% Clear temporary variables

clearvars filename delimiter formatSpec fileID dataArray ans;
set(data, 'OutputName', 'Stress');
set(data, 'OutputUnit', 'MPa');
set(data, 'InputName', 'Strain');
set(data, 'InputUnit', 'number');
set(data, 'Tstart',1, 'TimeUnit', 's');%0.0024


%% PART - II MODEL


%% Parameter definition

parameters = {1,1,1,1,1,1,1,1,1,1,1};

%% Order of the model % Order 
% Vector with three entries [Ny Nu Nx], specifying the number of model
% outputs Ny, the number of inputs Nu, and the number of states Nx

order         = [1 1 2];  
                           
                           
%% State variable initialization


rho_m0 = 1e11;
rho_f0 = 1e10; 

initial_states = [rho_m0; rho_f0];

%% Model definition


TIVmodel = idnlgrey('TIVSumDescrete', order, parameters, initial_states);   % for continuous data

TIVmodel.Algorithm.MaxIter = 100;


%%  Fixiing of parameters             

setpar(TIVmodel, 'Fixed', {false false false false false true true true true true true});

%% Fixing of initial state variables

setinit(TIVmodel,'Fixed',{true true});

%% To fix algorithm options. If we dont fix then it uses auto mode to select the algorithm options

 set(TIVmodel,'TimeUnit','s');
 set(TIVmodel,'CovarianceMatrix','Estimate');
 TIVmodel.Algorithm.SearchMethod = 'auto'; %'lm',gna,gn,grad,auto
 
               
TIVmodel.Algorithm.SimulationOptions.Solver = 'auto';


%% Setting parameters' values and limits


TIVmodel.Parameters(1).Value = 7e7;             %k1
TIVmodel.Parameters(2).Value = 10;              %k2
TIVmodel.Parameters(3).Value = 70;              %k3
TIVmodel.Parameters(4).Value = 81;             %kL
TIVmodel.Parameters(5).Value= 15e-6;             %Ls
TIVmodel.Parameters(6).Value = 25e-6;           %L0
TIVmodel.Parameters(7).Value = 2.96;            %M
TIVmodel.Parameters(8).Value = 2.86e-10;        %b
TIVmodel.Parameters(9).Value = 1/3;             %alpha
TIVmodel.Parameters(10).Value = 89; % 93.86; 262         %sigma_i
TIVmodel.Parameters(11).Value = 27e3;           %G MPa


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

                        
