% Script for computing p-safety for PDMP
% Dependencies: 1) Yalmip 2) MOSEK
% Authors: Rahul Misra, Sina Hassani, Rafal Wisniewski, and Manuela L.
% Bujorianu

clc
close all
clear 
yalmip('clear')

% Sinas IEEE 39 bus

% % % % % % % % % % % % % % % % % % % % Data of  Area 1
B1=0.3417;
R1=3;
Tg1=0.08;
Tt1=0.4;
H1=0.1677/2;
D1=0.084;
T12=0.2;
T13=0.25;
Kp1=0.33;
Ki1=-0.37;
% % % % % % % % % % % % % % % % % % % % Data of  Area 2
B2=0.3823;
R2=2.73;
Tg2=0.06;
Tt2=0.44;
H2=0.2017/2;
D2=0.016;
T21=0.2;
T23=0.12;
Kp2=0.36;
Ki2=-0.34;
% % % % % % % % % % % % % % % % % % % % Data of  Area 3
B3=0.3696;
R3=2.82;
Tg3=0.07;
Tt3=0.3;
H3=0.1247/2;
D3=0.015;
T31=0.25;
T32=0.12;
Kp3=0.25;
Ki3=-1;
Sbase=100;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
A11=[-1/Tg1 0 -1/(R1*Tg1) 0;
    1/Tt1 -1/Tt1 0 0;
    0 1/(2*H1) -D1/(2*H1) -1/(2*H1);
    0 0 2*pi*(T12+T13) 0];
A22=[-1/Tg2 0 -1/(R2*Tg2) 0;
    1/Tt2 -1/Tt2 0 0;
    0 1/(2*H2) -D2/(2*H2) -1/(2*H2);
    0 0 2*pi*(T21+T23) 0];
A33=[-1/Tg3 0 -1/(R3*Tg3);
    1/Tt3 -1/Tt3 0;
    0 1/(2*H3) -D3/(2*H3)];
A12=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 -2*pi*T12 0];

A13=[0 0 0;0 0 0;0 0 0;0 0 -2*pi*T13];
A21=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 -2*pi*T21 0];
A23=[0 0 0;0 0 0;0 0 0;0 0 -2*pi*T23];
A31=[0 0 0 0;0 0 0 0;0 0 0 1/(2*H3)];
A32=[0 0 0 0;0 0 0 0;0 0 0 1/(2*H3)];

b1=[1/Tg1;0;0;0];
b2=[1/Tg2;0;0;0];
b3=[1/Tg3;0;0];
c1=[0 0 B1 1];
c2=[0 0 B2 1];
% c3=[0 0 B3 1];
% c1=[0 0 1 0];
% c2=[0 0 1 0];
% c3=[0 0 1 0];
e1=[0;0;-1/(2*H1);0];
e2=[0;0;-1/(2*H2);0];
e3=[0;0;-1/(2*H3)]; %% positive or negative

A_p=[A11 A12 A13;
    A21 A22 A23;
    A31 A32 A33];
B_p=[b1 zeros(4,1) zeros(4,1);
    zeros(4,1) b2 zeros(4,1);
    zeros(3,1) zeros(3,1) b3];
C_p=[c1 zeros(1,4) zeros(1,3);
    zeros(1,4) c2 zeros(1,3);
    0 0 0 -1 0 0 0 -1 0 0 B3];
E_p=[e1 zeros(4,1) zeros(4,1);
    zeros(4,1) e2 zeros(4,1);
    zeros(3,1) zeros(3,1) e3];

% Now Controllers
A_c=zeros(3,3)-0*eye(3);
B_c=diag([Ki1,Ki2,Ki3]*0.2);
C_c=diag([1 1 1]);
D_c=diag([Kp1,Kp2,Kp3]*0.3);


% Augmenting the system
A_11=[A_p+B_p*D_c*C_p B_p*C_c;B_c*C_p A_c];
A_12=[B_p*D_c;B_c];
A_21=-[C_p*(A_p+B_p*D_c*C_p) C_p*(B_p*C_c)];
A_22=-C_p*(B_p*D_c);
A=[A_11 A_12;A_21 A_22];
A_13=[E_p;zeros(3,3)];
A_23=[C_p*E_p];
 C=zeros(3,14);
 C(1,3)=1;
 C(2,7)=1;
 C(3,11)=1;
D=zeros(3,3);
 E=[E_p;zeros(3,3);C_p*E_p];

% Overall dynamics with x vectr [x, e, \tau, s] 
A = blkdiag(A,1, zeros(3,3));

%% 
x = sdpvar(size(A,1),1); % the states
sdpvar p 
D = 2; % the degree of all involved SOS polynomials

% Linear Dynamics
F = A*x;

% Reset Map
R1 = [x(1:17); 0; x(15:17)]; % R1 = [x, e, 0, e]
R2 = [x(1:14); x(15:17) - x(19:21); x(18); zeros(3,1)]; % R1 = [x, e - s, \tau, zeros(3,1)]

% Max time in seconds
T_max = 1e-1;
tau = x(18);

%% We use 9 unknown (sum of squares) polynomials
[s1,coefs1] = polynomial(x,D);
[s2,coefs2] = polynomial(x,D);
[s3,coefs3] = polynomial(x,D);
[s4,coefs4] = polynomial(x,D);
[s5,coefs5] = polynomial(x,D);
[s6,coefs6] = polynomial(x,D);
[s7,coefs7] = polynomial(x,D);
[s8,coefs8] = polynomial(x,D);
[s9,coefs9] = polynomial(x,D); % SOS multiplier for boundary condition
[s10,coefs10] = polynomial(x,D); % SOS multiplier for boundary condition

% Choose your lambda
lambda = 1e0;

% The definition of three S, U, and E
% The state space  S = [S(i,x) \geq 0, i = 1,2] 
S1 = x(3) - 0.8;
S2 = x(7) - 0.8;
S3 = x(11) - 0.8;
S4 = -x(3) - 0.8;
S5 = -x(7) - 0.8;
S6 = -x(11) - 0.8;

% The initial set A = [A(i,x) \geq 0, i = 1,2]
A1 = x(3) - 0.8;
A2 = x(7) - 0.8;
A3 = x(11) - 0.8;
A4 = -x(3) - 0.2;
A5 = -x(7) - 0.2;
A6 = -x(11) - 0.2;

A7 = x(3) + 0.2;
A8 = x(7) + 0.2;
A9 = x(11) + 0.2;
A10 = -x(3) - 0.8;
A11 = -x(7) - 0.8;
A12 = -x(11) - 0.8;

% The forbidden state  U = [U(i,x) \geq 0, i = 1,2]
U1 = x(3) + 0.81;
U2 = x(7) + 0.81;
U3 = x(11) + 0.81;
U4 = -x(3) + 0.81;
U5 = -x(7) + 0.81;
U6 = -x(11) + 0.81;

% The two test functions for Subsystem#1 and Subsystem#2 h = (h_1,h_2)
[h1,coefh1,monh1] = polynomial(x,D);
[h2,coefh2,monh2] = polynomial(x,D);

% h1(R1) and h2(R2)
h1R1 = 0;
for iter = 1:length(monh1)
    new_mono = replace(monh1(iter), x, R1);   % or any substitution
    h1R1 = h1R1 + coefh1(iter)*new_mono;
end

h2R2 = 0;
for iter = 1:length(monh1)
    new_mono1 = replace(monh1(iter), x, R2);   % or any substitution
    h2R2 = h2R2 + coefh2(iter)*new_mono1;
end

% The infinitesmal generator consists of two components
Lh1 = jacobian(h1,x)*F + lambda*(h1R1-h1);
Lh2 = jacobian(h2,x)*F + lambda*(h2R2-h2);

% Constraints for the optimisation
% We use eight sos: s1, s2, s3, s4, s5, s6, s7, s8
constr = [sos(s1);sos(s2);sos(s3);sos(s4);sos(s5);sos(s6);sos(s7);sos(s8)];
% h >= 0 on the state space S 
constr = [constr; sos(h1 - s1*S1)];
constr = [constr; sos(h2 - s2*S2)];

% The infnitesmal generator -Lh >= 0 on the state space S 
constr = [constr; sos(-Lh1 - s3*S1)];
constr = [constr; sos(-Lh2 - s4*S2)];

% p - h >= 0 on A
constr = [constr; sos(p-h1 - s5*A1)];
constr = [constr; sos(p-h2 - s6*A2)];

% h - 1 >= 0 on U
constr = [constr; sos(h1-1 - s7*U1)];
constr = [constr; sos(h2-1 - s8*U2)];

% Boundary condition: h(R1) <= h(R2)
G1 = tau - T_max; % tau is the timer state
G2 = -tau + T_max;
constr = [constr; sos(+h2R2 - h1R1 - s9*G1 - s10*G2)];

% p >= 0
constr = [constr; p>=0; p<=1];

% % h1(0) == 0 and h2(0) == 0
% constr = [constr; coefh1(1) == 0; coefh2(1) == 0];

% SOLUTION
params = [coefs1;coefs2;coefs3;coefs4;coefs5;coefs6;coefs7;coefs8;coefh1;coefh2;p];
params = [params; coefs9];
options = sdpsettings('solver','mosek');
[sol, v, Q] = solvesos(constr,p,options,params);

% FeasibilityConstraints = [constr, p <= 0.9999*value(p)];
% [sol, v, Q] = solvesos(FeasibilityConstraints,[],options,params);
% 
% value(p)      
min([eig(Q{1}),eig(Q{2}),eig(Q{3}),eig(Q{4}),eig(Q{5}),eig(Q{6}),eig(Q{7}),eig(Q{8})])
min([eig(Q{9}),eig(Q{10}),eig(Q{11}),eig(Q{12}),eig(Q{13}),eig(Q{14}),eig(Q{15}),eig(Q{16})]) 

Br1 = clean(replace(h1,coefh1,double(coefh1)),1e-8);
sdisplay(Br1)

Br2 = clean(replace(h2,coefh2,double(coefh2)),1e-8);
sdisplay(Br2)
