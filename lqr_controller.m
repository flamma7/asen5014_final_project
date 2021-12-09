close all; clear all; clc;

n = sqrt(398600 / 6778^3 );

A = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];
B = zeros(6,3);
B(4:end,:) = eye(3);
C = zeros(3,6);
C(:, 1:3) = eye(3);
D = zeros(3);

P = ctrb(A,B);
rank(P); % System is fully controllable

% Check controllability
O = obsv(A,C);
rank(O); % System is fully observable

OLsys = ss(A,B,C,D); 

%% Step 2. Define reference input profiles and constraints 
umax = 1;
tvec = 0:0.01:60;
rhistvec = zeros(3, length(tvec));
rhistvec(1,:) = sign(double(tvec > 1 & tvec < 30)); % Set reference input to 1 from t=1 to t=30

%% Set up Integral control
Aaug = [A zeros(6,3); -C zeros(3,3)]; 
Baug = [B; 
        zeros(size(C,1),size(B,2))];
Faug = [zeros(size(B)); 
          eye(3);
         zeros(6,3)];
Caug = [C, zeros(3,3)];
Daug = zeros(size(Caug,1),size(Baug,2));

%%Assess reachability: 
rank(ctrb(Aaug,Baug)); %should be = 9

% Set poles for K and L seperately
% despoles_K = -[1 2 3 4 5 6 7 8 9]*0.2;
despoles_L = -[1 2 3 4 5 6]*0.2;

% Q is nxn
% R is mxm
% Strategy: define the diagonal vector I want then use diag() command
ai = 1/9;
bi = 1/3;
pos_max = 5;
vel_max = 5;
e_max = 5;
u_max = 1;
rho = 0.05;
q_diag = [ai^2/pos_max, ai^2/pos_max, ai^2/pos_max, ai^2/vel_max, ai^2/vel_max, ai^2/vel_max, ai^2/e_max, ai^2/e_max, ai^2/e_max];
q_diag(3) = 2*q_diag(3);
r_diag = [bi^2 / u_max, bi^2/u_max, bi^2/u_max];
Q = diag(q_diag);
R = rho * diag(r_diag);
Rtemp = diag(r_diag);

AugOLsys = ss(Aaug,Baug,Caug,Daug); 

[Kaug,Waug,clEvalsAug] = lqr(AugOLsys,Q,R);

% PLOT THE ROOT LOCUS
figure(), hold on
set(gca,'Color','k')
%%plot open loop poles
OLpoles = eig(A);
plot(OLpoles,'gx','LineWidth',2.5,'MarkerSize',9)
rhovals = logspace(-6,3,1000);
CLEvalsAug = nan(length(q_diag),length(rhovals));
for rr=1:length(rhovals)
    [~,~,CLEvalsAug(:,rr)] = lqr(AugOLsys,Q,rhovals(rr)*Rtemp);
    plot(real(CLEvalsAug(:,rr)),imag(CLEvalsAug(:,rr)),'w.')
end
title("Root Locus of Integral Error LQR Closed Loop Poles");

clEvalsAug

% Kaug = place(Aaug,Baug,despoles_K); 
L=(place(A.',C.', despoles_L)).';

% Add observer
AaugCLO = [(Aaug - Baug*Kaug) Baug*Kaug(:,1:6);
    zeros(6,9) (A-L*C)];
BaugCLO = Faug;
CaugCLO = [C zeros(3,9)];
DaugCLO = zeros(size(CaugCLO,1),size(BaugCLO,2));

CLaugsys2 = ss(AaugCLO,BaugCLO,CaugCLO,DaugCLO); 

%XCLO_IC = 0*ones(15,1); %zero initial error
XCLO_IC = zeros(15,1);
%XCLO_IC(10:15,1) = 0.1; %non-zero initial error

%% Step 5. Check that closed-loop system specs met; change despoles o'wise
%%get response to first reference input profile:
[Y_CL1,~,X_CL] = lsim(CLaugsys2,rhistvec,tvec,XCLO_IC);
U_CL = -[Kaug, Kaug(:,1:6)]*X_CL';
X_CL = X_CL';

%% PLOT ACTUATOR EFFORTS
figure() 
subplot(131)
plot(tvec, U_CL(1,:),'r') 
ax = gca;
ax.FontSize = 16; 
hold on
plot(tvec,umax*ones(size(tvec)),'k--') 
plot(tvec,-umax*ones(size(tvec)),'k--') 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('Thrust (N)', 'FontSize', 24) 
title('x thruster vs time', 'FontSize', 24) 

subplot(132)
plot(tvec, U_CL(2,:),'r') 
ax = gca;
ax.FontSize = 16; 
hold on
plot(tvec,umax*ones(size(tvec)),'k--') 
plot(tvec,-umax*ones(size(tvec)),'k--') 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('Thrust (N)', 'FontSize', 24) 
title('y thruster vs time', 'FontSize', 24) 

subplot(133)
plot(tvec, U_CL(3,:),'r') 
ax = gca;
ax.FontSize = 16; 
hold on
plot(tvec,umax*ones(size(tvec)),'k--') 
plot(tvec,-umax*ones(size(tvec)),'k--') 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('Thrust (N)', 'FontSize', 24) 
title('z thruster vs time', 'FontSize', 24) 

%% PLOT STATES

figure()
% X plots
subplot(321), hold on
plot(tvec, X_CL(1,:),'r') 
xline(1)
xline(11)
xline(30)
xline(40)
ax = gca;
ax.FontSize = 16; 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('x [km]', 'FontSize', 24) 
title('x (radial) vs time', 'FontSize', 24) 

subplot(322), hold on
plot(tvec, X_CL(4,:),'r') 
xline(1)
xline(11)
xline(30)
xline(40)
ax = gca;
ax.FontSize = 16; 
xlabel('t (secs)') 
ylabel('$\dot{x}$ [km/s]', 'Interpreter', 'latex', 'FontSize', 24) 
title('$\dot{x}$ (radial) vs time', 'Interpreter', 'latex', 'FontSize', 24)
legend('States','10 Second Timing Objective', 'FontSize', 18)


% Y plots
subplot(323), hold on
plot(tvec, X_CL(2,:),'r') 
xline(1)
xline(11)
xline(30)
xline(40)
ax = gca;
ax.FontSize = 16; 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('y [km]', 'FontSize', 24) 
title('y (along-track) vs time', 'FontSize', 24)

subplot(324), hold on
plot(tvec, X_CL(5,:),'r') 
xline(1)
xline(11)
xline(30)
xline(40)
ax = gca;
ax.FontSize = 16; 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('$\dot{y}$ [km/s]', 'Interpreter', 'latex', 'FontSize', 24) 
title('$\dot{y}$ (along-track) vs time', 'Interpreter', 'latex', 'FontSize', 24)


% Z plots
subplot(325), hold on
plot(tvec, X_CL(3,:),'r') 
xline(1)
xline(11)
xline(30)
xline(40)
ax = gca;
ax.FontSize = 16; 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('z [km]', 'FontSize', 24) 
title('z (cross-track) vs time', 'FontSize', 24)

subplot(326), hold on
plot(tvec, X_CL(6,:),'r') 
xline(1)
xline(11)
xline(30)
xline(40)
ax = gca;
ax.FontSize = 16; 
xlabel('t (secs)', 'FontSize', 24) 
ylabel('$\dot{z}$ [km/s]', 'Interpreter', 'latex', 'FontSize', 24) 
title('$\dot{z}$ (cross-track) vs time', 'Interpreter', 'latex', 'FontSize', 24)