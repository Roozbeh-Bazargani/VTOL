%% Linearization

%state variables
syms y yd p pd r rd thld thrd X Xd Y Yd Z Zd
%input variables
syms Vl Vr Vb phil phir

%constants
% ydd
l3 = 0.5;
Jy = 0.25934;
Kdrear = 0;
Kdprop = 0;
Kprop = 2.5*9.81/10;
Fl = Kprop * Vl;
Fr = Kprop * Vr;
% pdd
l1 = 0.05;
l2 = 0.8;
Jp = 0.094625;
Krear = 0.015;
Fb = Krear * Vb;
% rdd
Jr = 0.197212;
KD = 0.004;
Kt = 0.1901;
Ra = 0.29;
Ke = Kt;
% thmdd
bm = 0;
Jprop = 2.55e-1;
% zdd
Mp = 4;
g = 9.8;
W = Mp * g;

A_nonlinear = [
    yd;
    0;
    pd;
    0;
    rd;
    0;
    -(bm + (Kt*Ke)/Ra)/Jprop*thld - KD/Jprop*thld^2;
    -(bm + (Kt*Ke)/Ra)/Jprop*thrd - KD/Jprop*thrd^2;
    Xd;
    0;
    Yd;
    0;
    Zd;
    0;
];

            
B_nonlinear = [
    0; 
    (l3/Jy)*(Fl*sin(phil) - Fr*sin(phir)) + (Kdprop/Jy)*(Vl - Vr) + (Kdrear/Jy)*Vb; %l3 - Jy - *Fl - *Fr - Kdprop - Kdrear
    0;
    (l1/Jp)*(Fl*cos(phil) + Fr*cos(phir)) - (l2/Jp)*Fb; %l1 - Jp - *Fl - *Fr - l2 - *Fb
    0;
    (l3/Jr)*(Fl*cos(phil) - Fr*cos(phir));
    (Kt/(Ra*Jprop))*Vl;
    (Kt/(Ra*Jprop))*Vr;
    0;
    (Fl*sin(phil) + Fr*sin(phir)) / Mp;
    0;
    0;
    0;
    Fl * cos(phil) + Fr * cos(phir) + Fb - W;
];
       
states = [y yd p pd r rd thld thrd X Xd Y Yd Z Zd];
inputs = [Vl Vr Vb phil phir];
%Sop = solve(transpose(A_nonlinear + B_nonlinear), [states inputs]);
%op = vpa([S.y S.yd S.p S.pd S.r S.rd S.thld S.thrd S.X S.Xd S.Y S.Yd S.Z S.Zd S.Vl S.Vr S.Vb S.phil S.phir]',5);
op = [ 0, 0, 0, 0, 0, 0, 22.833, 22.833, 0, 0, 0, 0, 0, 0, 7.5217, 7.5217, 153.73, 0, 0];
states_op = op(1:14);
inputs_op = op(15:end);

As = jacobian(A_nonlinear, states);
A = double(vpa(simplify(subs(As, states, states_op))));

Bs = jacobian(B_nonlinear, inputs);
B = double(vpa(simplify(subs(Bs, inputs, inputs_op))));
% sensors
C = [1 zeros(1, 13);
    0 0 1 zeros(1, 11);
    0 0 0 0 1 zeros(1, 9);
    zeros(1, 8) 1 0 0 0 0 0;
    zeros(1, 12) 1 0];

%% clearing 
save('space state', 'A', 'B', 'C')
clear;
load('space state')


%% state transformation matrix
syms t;
phi = vpa(expm(A * t), 3);
% y yd p pd r rd thld thrd X Xd Y Yd Z Zd
initial_state = [0, 0, 0, 0, 0, 0, 220.833, 220.833, 500, 0, 0, 0, 0, 0,];
%initial_state_response_linear = vpa(transpose(phi * transpose(initial_state)), 3)
%initial_state_response_nonlinear = vpa(subs(transpose(A_nonlinear), states, initial_state), 3)


%% controllabilty and observablity
syms s
phiC = ctrb(A, B);
phiO = obsv(A, C);
phiC_rank = rank(ctrb(A, B))
PhiO_rank = rank(obsv(A, C))
  
Tc = [phiC(1:end, 1:12), [zeros(10, 1); 1; 0;0;0], [zeros(11, 1); 1; 0;0]];
A_hat = inv(Tc) * A * Tc;
B_hat = inv(Tc) * B;
C_hat = C * Tc;
A_hat_c = A_hat(1:12, 1:12);
B_hat_c = B_hat(1:12, 1:end);
C_hat_c = C_hat(1:end, 1:12);


C_o = obsv(A_hat_c, C_hat_c);
Tc_o = [C_o(1:10, 1:12); [zeros(1, 10), 1, 0]; [zeros(1, 10), 0, 1]];
A_hat_c_hat = Tc_o * A_hat_c * inv(Tc_o);
B_hat_c_hat = Tc_o * B_hat_c;
C_hat_c_hat = C_hat_c * inv(Tc_o);
A_hat_c_o = A_hat_c_hat(1:10, 1:10);
B_hat_c_o = B_hat_c_hat(1:10, 1:end);
C_hat_c_o = C_hat_c_hat(1:end, 1:10);
A_hat_c_no = A_hat_c_hat(11:12, 11:12);
B_hat_c_no = B_hat_c_hat(11:12, 1:end);
C_hat_c_no = C_hat_c_hat(1:end, 11:12);


A_hat_nc = A_hat(13:14, 13:14);
B_hat_nc = B_hat(13:14, 1:end);
C_hat_nc = C_hat(1:end, 13:14);
nC_o = obsv(A_hat_nc, C_hat_nc);
Tnc_o = [1, 0; 0, 1];
A_hat_nc_no = Tnc_o * A_hat_nc * inv(Tnc_o);
B_hat_nc_no = Tnc_o * B_hat_nc;
C_hat_nc_no = C_hat_nc * inv(Tnc_o);


%% full observer
scale = -1;
%A_full_obsv = [A(1:6, 1:6), A(1:6, 9:10), A(1:6, 13:14);A(9:10, 1:6), A(9:10, 9:10), A(9:10, 13:14);
    %A(13:14, 1:6), A(13:14,9:10), A(13:14, 13:14)]
    A_full_obsv = [A(1:10, 1:10), A(1:10, 13:14); A(13:14, 1:10), A(13:14, 13:14)];
% A_full_obsv = [A(1:6, 1:6), zeros(6, 2),A(1:6, 9:10), zeros(6, 2),A(1:6, 13:14); zeros(2, 14);A(9:10, 1:6), zeros(2, 2),  A(9:10, 9:10), zeros(2, 2) , A(9:10, 13:14);
%     zeros(2, 14); A(13:14, 1:6), zeros(2,2), A(13:14,9:10), zeros(2,2),A(13:14, 13:14)]
%C_full_obsv = [C(1:end, 1:6), C(1:end, 9:10), C(1:end, 13:14)]
C_ro = [1 zeros(1, 13);
    0 0 1 zeros(1, 11);
    0 0 0 0 1 zeros(1, 9);
    zeros(1, 8) 1 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0;
    zeros(1, 12) 1 0];
C_ro1 = [C_ro(1:end, 1:10), C_ro(1:end, 13:14)];
% L = transpose(acker(transpose(A), transpose(C))
% 7, 8, 11, 12
% thetadl, thetadr, Y, Yd
L = transpose(place(transpose(A_full_obsv), transpose(C_ro1), [-10, -11, -12, -13, -14, -15, -16, -17, -18, -19 -20 -21]));
L = [L(1:10, 1:end); zeros(2,7); L(11:end, 1:end)];
Ac = [A(1:10,1:10), A(1:10, 13:end);A(13:end,1:10), A(13:end, 13:end)];
Bc = [B(1:10,1:end);B(13:end,1:end)];
Cc = [1 zeros(1, 11);
    0 0 1 zeros(1, 9);
    0 0 0 0 1 zeros(1, 7);
    zeros(1, 8) 1 0 0 0;
    zeros(1, 10) 1 0];
%rank(ctrb(Ac, Bc))
pK = [-10 -11 -12 -13 -14 -15 -16 -17 -18 -19 -20 -21] / 10;
Ktmp = place(Ac, Bc, pK);
K = [Ktmp(1:end, 1:10), zeros(5, 2), Ktmp(1:end, 11:end)];

%% reduced state observer
C_ro = [1 zeros(1, 13);
    0 0 1 zeros(1, 11);
    0 0 0 0 1 zeros(1, 9);
    zeros(1, 8) 1 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0;
    zeros(1, 12) 1 0];
C_ro = [C_ro(1:end, 1:10), C_ro(1:end, 13:14)];
l = 7;
n = 12;
D = [-10, zeros(1, 4); 0, -11, zeros(1, 3); zeros(1, 2), -12, zeros(1, 2); zeros(1, 3), -13, zeros(1, 1); zeros(1, 4), -14];
cq = [zeros(l, n - l), eye(l)]; %C * Q
B_ro = [B(1:10, 1:end); B(13:14, 1:end)];
% Q = [zeros(1, 5), 1, zeros(1, 4); 1, zeros(1,9); zeros(1, 6), 1, zeros(1, 3);
%     0, 1, zeros(1, 8); zeros(1, 7), 1, zeros(1, 2); 0, 0, 1, zeros(1, 7);
%     zeros(1, 8), 1, 0; zeros(1, 3), 1, zeros(1, 6); zeros(1, 9), 1; zeros(1, 4), 1, zeros(1, 5)];
% 
% Atil = inv(Q) * A_full_obsv * Q
% 
% Atil_11 = Atil(1:5, 1:5)
% Atil_12 = Atil(1:9, 10:10)
% Atil_21 = Atil(10:10, 1:9)
% Atil_22 = Atil(10, 10)
% 
% Ltil = (D - Atil_11) * pinv(Atil_21)
% 
Q = [zeros(1, 5), 1, zeros(1, 6); zeros(1, 1), 1, zeros(1, 10); zeros(1, 6), 1, zeros(1, 5); 1, zeros(1, 11); zeros(1, 7), 1, zeros(1, 4); zeros(1, 3), 1, zeros(1, 8); zeros(1, 9), 1, zeros(1, 2); zeros(1, 10), 1, zeros(1, 1); zeros(1, 8), 1, zeros(1, 3); zeros(1, 4), 1, zeros(1, 7); zeros(1, 11), 1; zeros(1, 2), 1, zeros(1, 9)];
A_reduced_obsv = [A(1:10, 1:10), A(1:10, 13:14); A(13:14, 1:10), A(13:14, 13:14)];
Atil = inv(Q) * A_reduced_obsv * Q;

Atil_11 = Atil(1:5, 1:5);
Atil_12 = Atil(1:5, 6:end);
Atil_21 = Atil(6:end, 1:5);
Atil_22 = Atil(6:end, 6:end);

Ltil = (D - Atil_11) * pinv(Atil_21);
L = [eye(n - l), Ltil] * inv(Q);
FF = inv([C_ro; L]);
F1 = FF(1:12, 1:7);
F2 = FF(1:12, 8:end);
R = L * B_ro;

T = Atil_12 + Ltil * Atil_22 - D*Ltil;

Kro = [K(1:end,1:10) K(1:end,13:end)];



%% Controller state feedback
Ac = [A(1:10,1:10), A(1:10, 13:end);A(13:end,1:10), A(13:end, 13:end)];
Bc = [B(1:10,1:end);B(13:end,1:end)];
Cc = [1 zeros(1, 11);
    0 0 1 zeros(1, 9);
    0 0 0 0 1 zeros(1, 7);
    zeros(1, 8) 1 0 0 0;
    zeros(1, 10) 1 0];
%rank(ctrb(Ac, Bc))
pK = [-10 -11 -12 -13 -14 -15 -16 -17 -18 -19 -20 -21] / 10;
Ktmp = place(Ac, Bc, pK);
K = [Ktmp(1:end, 1:10), zeros(5, 2), Ktmp(1:end, 11:end)];

%% Static pre-filter
Gcl = -[0 0 0 0 1] * Cc * inv(Ac-Bc*Ktmp) * Bc;
pGcl = pinv(Gcl);

%% Dynamic pre-filter
Cci = [zeros(1, 10) 1 0];
Anew = [Ac, zeros(12, 1); -Cci, zeros(1, 1)];
Bnew = [Bc; zeros(1, 5)];
Cnew = [Cci, zeros(1, 1)];
pKi = [-10 -11 -12 -13 -14 -15 -16 -17 -18 -19 -20 -21, -22] / 1;
alpha = poly(pKi);
alphaA = polyval(alpha, Anew);
k3tmp = pinv(ctrb(Anew, Bnew)) * alphaA;
k3tmp = k3tmp(61:end, 1:end);
K3 = [k3tmp(1:end, 1:10), zeros(5, 2), k3tmp(1:end, 11:end)];
% k3tmp = place(Anew, Bnew, pKi);
%K3 = [k3tmp(1:end, 1:10), zeros(5, 2), k3tmp(1:end, 11:end)]

%% optimized Control
Ac = [A(1:10,1:10), A(1:10, 13:end);A(13:end,1:10), A(13:end, 13:end)];
Bc = [B(1:10,1:end);B(13:end,1:end)];
Q=diag([10,1,10,1,10,1,20,20,10,1,100,1]);
R=diag([10,10,10,1,1]);
[X,L,G] = care(Ac,Bc,Q,R,zeros(12,5),eye(12));
K = [G(1:end,1:10),zeros(5,2),G(1:end,11:end)];

 







