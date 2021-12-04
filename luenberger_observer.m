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
rank(P) % System is fully controllable

% Check controllability
O = obsv(A,C);
rank(O) % System is fully observable

OLsys = ss(A,B,C,D); 

%% Step 2. Define reference input profiles and constraints 
umax = 1;
tvec = 0:0.01:60;
rhistvec = zeros(3, length(tvec));
rhistvec(1,:) = sign(double(tvec > 1 & tvec < 30)); % Set reference input to 1 from t=1 to t=30

%% Set up Luenberger observer
Aaug = [A zeros(6,3); -C zeros(3,3)]; 
Baug = [B; 
        zeros(size(C,1),size(B,2))];
Faug = [zeros(size(B)); 
          eye(3);
         zeros(6,3)];
Caug = [C, zeros(3,3)];
Daug = zeros(size(Caug,1),size(Baug,2));

%%Assess reachability: 
rank(ctrb(Aaug,Baug)) %should be = 9

% Set poles for K and L seperately
despoles_K = -[1 2 3 4 5 6 7 8 9]*0.2;
despoles_L = -[1 2 3 4 5 6]*0.2;

Kaug = place(Aaug,Baug,despoles_K); 
L=(place(A.',C.', despoles_L)).';

AaugCLO = [(Aaug - Baug*Kaug) Baug*Kaug(:,1:6);
    zeros(6,9) (A-L*C)];
BaugCLO = Faug;
CaugCLO = [C zeros(3,9)];
DaugCLO = zeros(size(CaugCLO,1),size(BaugCLO,2));

CLaugsys2 = ss(AaugCLO,BaugCLO,CaugCLO,DaugCLO); 

%XCLO_IC = 0*ones(15,1); %zero initial error
XCLO_IC = 0.1*ones(15,1); %non-zero initial error

%% Step 5. Check that closed-loop system specs met; change despoles o'wise
%%get response to first reference input profile:
[Y_CL1,~,X_CL] = lsim(CLaugsys2,rhistvec,tvec,XCLO_IC);
U_CL = -[Kaug, Kaug(:,1:6)]*X_CL';
X_CL = X_CL';

%% PLOT ACTUATOR EFFORTS
figure() 
subplot(131)
plot(tvec, U_CL(1,:),'r') 
hold on
plot(tvec,umax*ones(size(tvec)),'k--') 
plot(tvec,-umax*ones(size(tvec)),'k--') 
xlabel('t (secs)') 
% ylabel('') 
title('x thruster vs time') 

subplot(132)
plot(tvec, U_CL(2,:),'r') 
hold on
plot(tvec,umax*ones(size(tvec)),'k--') 
plot(tvec,-umax*ones(size(tvec)),'k--') 
xlabel('t (secs)') 
% ylabel('') 
title('y thruster vs time') 

subplot(133)
plot(tvec, U_CL(3,:),'r') 
hold on
plot(tvec,umax*ones(size(tvec)),'k--') 
plot(tvec,-umax*ones(size(tvec)),'k--') 
xlabel('t (secs)') 
% ylabel('') 
title('z thruster vs time') 

%% PLOT STATES

figure()
% X plots
subplot(321), hold on
plot(tvec, X_CL(1,:),'r') 
xlabel('t (secs)') 
ylabel('x [km]') 
title('x (radial) vs time') 

subplot(322), hold on
plot(tvec, X_CL(4,:),'r') 
xlabel('t (secs)') 
ylabel('$\dot{x}$ [km/s]', 'Interpreter', 'latex') 
title('$\dot{x}$ (radial) vs time', 'Interpreter', 'latex')


% Y plots
subplot(323), hold on
plot(tvec, X_CL(2,:),'r') 
xlabel('t (secs)') 
ylabel('y [km]') 
title('y (along-track) vs time')

subplot(324), hold on
plot(tvec, X_CL(5,:),'r') 
xlabel('t (secs)') 
ylabel('$\dot{y}$ [km/s]', 'Interpreter', 'latex') 
title('$\dot{y}$ (along-track) vs time', 'Interpreter', 'latex')


% Z plots
subplot(325), hold on
plot(tvec, X_CL(3,:),'r') 
xlabel('t (secs)') 
ylabel('z [km]') 
title('z (cross-track) vs time')

subplot(326), hold on
plot(tvec, X_CL(6,:),'r') 
xlabel('t (secs)') 
ylabel('$\dot{z}$ [km/s]', 'Interpreter', 'latex') 
title('$\dot{z}$ (cross-track) vs time', 'Interpreter', 'latex')