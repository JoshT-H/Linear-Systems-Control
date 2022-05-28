clear all; close all; clc;

% Please complete Question 1 here.
L  = [20]; % Inductance
c  = [0.1]; % Capacitance
R  = [4]; % Resistance

% Derivation of the state space equation and system matrices: Note T =
% transpose 

% 1.) [I',I'']T = [0,1; (-1/(L*c)),(-R/L)]*[I,I']T + [0, 1/L]T*v' 
% 2.)  x' = A*x + B*u

Ac = [0,1; (-1/(L*c)),(-R/L)]; % Continuous time A matrix
Bc = [0;1/L]; % Continuous time B matrix

%Obtained Values
% Ac = [0,1;-0.5,-0.2];
% Bc = [0;0.05];

% Please complete Question 2 here.
dt = [2*10^-3]; % Sampling rate.
T  = [20]; % Simulation duration.

% Method 1
% Continous-time to discrete-time conversion using Euler Approximation [1]
%A = eye(2)+ dt*Ac  
%B = dt* Bc

%Obtained Values Mehtod 1
%B = [1,0.002;-1.0e-03,0.9996]
%A = [0;1.0e-04]

%Method 2
% Contionous-time t0 discrete-time conversion using Intergral Approximation [1]
A = expm(Ac*dt); % Discrete time A matrix
B = Ac\(A-eye(2))*Bc; %Discrete time B matrix.

% Obtained Values Method 2
% A = [0.1,0.002;-9.998e-04,0.999]
% B = [9.999e-08;9.998e-05]

u  = [1]; % Control input.
x  = [0;0]; % Initial state.
S  = 10000;  % Number of simulation steps.
X  = zeros(2,S); % Matrix for storing data

for s = 1:S
    x = A*x + B*u;
    X(:,s) = x;
end

% Please complete Question 3 here.
% The Observer matrix gives us the relationship between the input x and the
% output of the system i.e yk= Cxk:  Since xk = [I,I'], an Observer Matrix
% of C = [1,0] will give use the current output relation we desire. 

C = [1 0]; % Observer ('C') matrix
H = [C; C*A]; % Observability matrix

% Obtained H values 
%H = [1,0;1,0.002000000000000]

% To calculate the limit of the poles;
% 1.) The transfer function G is set zero.
% 2.) We expland the Right Hand Side of the expression
% 3.) We will then obtain a second order polynomial expression for z
% 4.) Apply the quadratic function z = (-b+/-(b^2-(4*a*c)^1/2)/2a to obtain
% the poles of z.

syms z % Create a variables z 
G = C*inv((z*eye(2))-A)*B;% Transfer Function of Discrete LTI Systems
poles = round(poles(G),4); % The pole function simpy applies the steps metioned above to solve for z 

% The sort function sort the element of poles in accending order. 
z = [sort(poles)]; % Vector containing poles in order of increasing size. 

% Obtained Values
% z = [0.9998 - 0.0014i, 0.9998 + 0.0014i]

% References 
% 1.) http://eceweb1.rutgers.edu/~gajic/solmanual/slides/chapter8_DIS.pdf

% Please do not change the template code after this line.
save(['rlc_',number,'_',firstname,'_',surname]);

