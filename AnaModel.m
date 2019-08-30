%% Analytical model of 3 cross 3 phase machine. 

%% Defining Variables
syms Ll % leakage inductance
syms Lms % mutual inductance
syms Llm % leakage mutual between two winding sets. 
syms rs % resistance
syms omg % 2*pi*fs where fs is the synchronous frequency
syms t % time
syms ('Va1(t)','Va2(t)','Va3(t)','Vb1(t)','Vb2(t)','Vb3(t)','Vc1(t)','Vc2(t)','Vc3(t)') % Phase voltage vectors for three winding phases Va1 = Vm sin(omgt), Vb1 = Vm sin(omgt-2pi/3)
syms ('ia1(t)','ia2(t)','ia3(t)','ib1(t)','ib2(t)','ib3(t)','ic1(t)','ic2(t)','ic3(t)') % currents. 
syms Phi_pm % permanent magnet flux

%% Defining self and mutual inductances
L1 = Ll + 1.5*Lms; 
L2 = (sqrt(3)/2)*Lms + Llm;
L3 = 1.5*Lms+2*Llm; 
%% Inductance Matrices
Mat1 = L1*eye(3);
Mat2 = [L2 -L2 0;
        0 L2 -L2;
        -L2 0 L2];
Mat3 = [0 -L3 0;
        0 0  -L3;
        -L3 0 0];
%% Stationary Frame of reference transformation matrices
T1 = 2/3*[ 1   -0.5       -0.5;
           0   -sqrt(3)/2  sqrt(3)/2;
           0.5  0.5        0.5];
T2 = 2/3*[sqrt(3)/2 -sqrt(3)/2 0;
          -0.5 -0.5 1;
          0.5 0.5 0.5];
T3 = 2/3*[0.5 -1 0.5;
         -sqrt(3)/2 0 sqrt(3)/2;
         0.5 0.5 0.5];
%% Alpha Beta Current Vectors. 
I1 = [ia1(t); ib1(t);0];
I2 = [ia2(t); ib2(t);0];
I3 = [ia3(t); ib3(t);0];

mech_tht= pi/6; 

Phi_1 = Mat1*inv(T1)*I1 +Mat2*inv(T2)*I2 +Mat3*inv(T3)*I3 + Phi_pm*[sin(omg*t);sin(omg*t-2*pi/3); sin(omg*t-4*pi/3)];
Phi_2 = Mat1 *inv(T2)*I2 + Mat2.'*inv(T1)*I1 +Mat2*inv(T3)*I3 + Phi_pm*[sin(omg*t-mech_tht);sin(omg*t-(2*pi/3+mech_tht)); sin(omg*t-(4*pi/3+mech_tht))];
Phi_3 = Mat1 *inv(T3)*I3 + Mat2.'*inv(T2)*I2 +Mat3.'*inv(T1)*I1 + Phi_pm*[sin(omg*t-mech_tht);sin(omg*t-(2*pi/3+mech_tht)); sin(omg*t-(4*pi/3 + mech_tht))];


V1 = rs*inv(T1)*I1 + diff(Phi_1);
V2 = rs*inv(T2)*I2 + diff(Phi_2);
V3 = rs*inv(T3)*I3 + diff(Phi_3);

%% Transformation into Alpha Beta co-ordinate system or the stationary rotating frame of reference. 
Val_bt1 = T1*V1;
Val_bt2 = T2*V2;
Val_bt3 = T3*V3;
