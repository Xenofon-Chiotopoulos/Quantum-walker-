% This version used with Fit_AqP_v4.m

global  program_version ...
        linelength linecenter ...
        ntimes theta ...
        Kl Kr F ...
        SummaryOutput_qubit SummaryOutput_probs SummaryOutput_meanZ
    
% file information
%-------------------------
program_version =           'Program: Quantum-walk-1D-v2.m';
SummaryOutput_qubit =       'Faux-ton-1D-qubit-v2-op1.out';
SummaryOutput_probs =       'Faux-ton-1D-probs-v2-op1.out';
SummaryOutput_meanZ =       'Faux-ton-1D-meanZ-v2-op1.out';



% Fixed quantities
%-------------------------
linelength = 99;          linecenter = (linelength+1)/2;        % linelength to be odd
ntimes = 50;              
detector_point = 51;
detector_on = 0;
detect_superposition = 1;
detector_everywhere = 0;
detector_everywhere_iteration = 30;
theta = pi/2;

x_r = 1;
x_c = 0;
y_r = 1;
y_c = 0;

normalize = 1/(abs((complex(x_r,x_c))*(complex(x_r,x_c)))+abs((complex(y_r,y_c))*(complex(y_r,y_c))))^0.5;

init_qubit(1,1) = complex(x_r,x_c)*normalize;           % the qubit at z=0 at t=0 a
init_qubit(2,1) = complex(y_r,y_c)*normalize;           % the qubit at z=0 at t=0 b
%init_qubit(1,1) = (1/sqrt(2))*complex(1,0);           % the qubit at z=0 at t=0 a
%init_qubit(2,1) = (1/sqrt(2))*complex(0,1);           % the qubit at z=0 at t=0 b
init_qubit;

%init_qubit = [complex(1/sqrt(2),0); complex(0,1/sqrt(2))];

pauli_y_eig_1 = (1/sqrt(2))*[complex(0,-1);1];  %testing qubit for bells theroem
pauli_y_eig_2 = (1/sqrt(2))*[1;complex(0,-1)]; 

init_qubit_y_create = pauli_y_eig_1*init_qubit(1,1)+init_qubit(2,1)*pauli_y_eig_2;
init_qubit_y(1,1) = [0,1]*init_qubit_y_create;
init_qubit_y(2,1) = [1,0]*init_qubit_y_create;

%testing new detector that seperates superposed qubit
%-------------------------

test_matrix_first = [1;0];
test_matrix_second = [0;1];

init_qubit_first(1,1) = [1,0]*test_matrix_first;
init_qubit_first(2,1) = [0,1]*test_matrix_first;
init_qubit_second(1,1) = [1,0]*test_matrix_second;
init_qubit_second(2,1) = [0,1]*test_matrix_second;

init_qubit_first;               %for debugging
init_qubit_second;

% the operational arrays
%-------------------------

Kl = [0,-1,0,0];
Kr = [0,1,0,0];
Hadamard = (1.0/sqrt(2.0))*[1, 1 ; 1, -1];

Ry = [cos(theta/2),-sin(theta/2);sin(theta/2),cos(theta/2)];
Rz = [complex(cos(theta/2),-sin(theta/2)), 0;0,complex(cos(theta/2),sin(theta/2))];
RyRz = Ry*Rz;

RyRzReal = 0.5*[sin(theta),-cos(theta);-cos(theta),-sin(theta)];
RyRzImag = -0.5*[cos(theta),sin(theta);-sin(theta),cos(theta)];


% select operator
%-------------------------
F = RyRz;







