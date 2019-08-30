%% Script to get frequency response and couplings between phases

%% Parameters
rs =1.98; 
Lms = 641e-6;
Ll = 344.5200e-6;
Llm = 38e-6;
Phi_f =0.03;
f = 100; 
P = 6; 
Tm = 1; 
J = 0.275e-3;
%% Initial calculations
L0 = Ll + 1.5*Lms; 
L1 =  (sqrt(3)* Llm + 1.5 *Lms);
L2 =  2*Llm + 1.5*Lms;
omg = 2*pi*f;
%% Matrices
Mat1 = [L0 0 L1 0 L2 0; 
        0 L0 0 L1 0 L2;
        L1 0 L0 0 L1 0;
        0 L1 0 L0 0 L1;
        L2 0 L1 0 L0 0;
        0 L2 0 L1 0 L0];
Mat2 = [-rs -omg*L0 0 -omg*L1 0 -omg*L2;
        omg*L0 -rs omg*L1 0 omg*L2 0;
        0 -omg*L1 -rs -omg*L0 0 -omg*L1;
        omg*L1 0 omg*L0 -rs omg*L1 0;
        0 -omg*L2 0 -omg*L1 -rs -omg*L0;
        omg*L2 0 omg*L1 0 omg*L0 -rs];
D = Phi_f*omg*[ 1;0;1;0;1;0];
A = inv(Mat1)*Mat2;
B = inv(Mat1);
C = [1 0 0 0 0 0]; 
D = zeros(1,1);
%% Transfer function
% B(:,i) signifies Iq1/Vqi for all i = 1,2,3; 
[b, a] = ss2tf(A,B(:,1),C,D);
 H = tf(b,a);
 figure(1) 
 bode(H)
 hold on 
 grid on
% Temp1 = Mat2*[Idq1; Idq2; Idq3];