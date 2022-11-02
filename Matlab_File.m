clear all
close all
clc

%% DATA

m1 = 4.7764;
m2 = 0.5155;
k1 = 27950;
c1 = 7.152;

%% STATE-SPACE REPRESENTATION (WITHOUT MATRIX A)

B = [  0   ;
      1/m1 ;
       0   ;
     -1/m2  ];

E = [  0   ;
     -1/m1 ;
       0   ;
       0  ];

C = [1 0 0 0];
D = 0;
F = 0;

%% L2-GAIN PLOT

% Ranges of variation for c2 and k2:
k2=1500:1:3500;
c2=1:0.2:2;

for j=1:numel(c2)
    for i=1:numel(k2)
% Matrix A for different values of c2 and k2:
A(:,:,i,j) = [    0               1          0        0   ;
   (-k1-k2(i))/m1  (-c1-c2(j))/m1   k2(i)/m1  c2(j)/m1 ;
         0               0          0        1   ;
      k2(i)/m2         c2(j)/m2   -k2(i)/m2  -c2(j)/m2];
sys = ss(A(:,:,i,j),E,C,F);
% Vector of the peaks of the transfer functions between w and z:
gamma_vec(i,j) = getPeakGain(sys);
    end
end

figure
hold on
for j=1:numel(c2)
    plot(k2,gamma_vec(:,j),'LineWidth',2)
end
title('L2-gain of the passive system')
xlabel('k2 (N/m)') 
ylabel('L2-gain values') 
legend('c2=1 kg/s','c2=1.2 kg/s','c2=1.4 kg/s','c2=1.6 kg/s','c2=1.8 kg/s','c2=2 kg/s','Location','bestoutside')
hold off

%% STATIC FULL-STATE FEEDBACK DESIGN

% Considered values of c2 and k2:
k2new = 1619;
c2new = 1.717;
K_bar = 5e3; % limit on the norm of matrix K

% Matrix A:
Anew = [    0                 1            0         0   ;
       (-k1-k2new)/m1  (-c1-c2new)/m1   k2new/m1  c2new/m1 ;
            0                 0            0         1   ;
         k2new/m2         c2new/m2    -k2new/m2  -c2new/m2];

% Plant Dimensions
[n,~] = size(Anew); %[n x n]
[~,p] = size(E); %[n x p]
[m,~] = size(C); %[m x n]
[~,d] = size(F); %[m x d]

% -- Internal Stability: Eigenvalues Test
ev = real(eig(Anew));

disp('Internal Stability Test (open-loop)')
disp('*Open-loop eigenvalues real parts: ')
disp(['    ', mat2str(ev')]);
disp('*Eigenvalues test:')
if all(ev<0) %all eigenvalues have Negative real part
    disp(['    ' 'all eigenvalues have negative real part --> system is stable'])
else
    disp(['    ' 'at least one eigenvalue has positive real part --> system is unstable'])
end

% -- Controllability Test
ctr = ctrb(Anew,B);

disp('--------------------------------------------------------------------------')
disp('Controllability Test');
disp('*Controllability matrix test:');
%%%%%%%%%%%%%%%%%%%
if (rank(ctr) == n) %
%%%%%%%%%%%%%%%%%%%
    disp(['    ' 'ctrb matrix is full rank --> system is controllable'])
else
    disp(['    ' 'ctrb matrix is rank deficient --> system is NOT controllable'])
    disp('*plant dimension:');
    disp(['    ', num2str(n)]);
    disp('*ctrb matrix rank:');
    disp(['    ', num2str(rank(ctr)), newline]);
end

% -- Feedback Design
disp('--------------------------------------------------------------------------')
disp('*STATIC FULL-STATE FEEDBACK DESIGN*')

% -- Clear the internal memory of YALMIP
yalmip('clear')

% -- Set the parameters of the LMI solver
opts = sdpsettings('solver', 'mosek', 'verbose', 0);

% -- Optimization variables

fake_zero = 1e-6;
W = sdpvar(n,n,'symmetric');
X = sdpvar(p,n);
gamma = sdpvar(1,1);
rho = sdpvar(1,1);

% -- Constraints

constr = [ rho >= fake_zero; ...
               W >= rho*eye(n); ...
               [ (Anew*W+B*X)+(Anew*W+B*X)',        E,        (C*W+D*X)'  ; ...
                             E',              -gamma*eye(d),      F'      ; ...
                         (C*W+D*X),                 F,      -gamma*eye(m) ] <= -fake_zero*eye(n+d+m); ...
               [K_bar*rho*eye(n),      X'      ; ...
                       X,      K_bar*rho*eye(p) ] >= 0 ];


% -- Solve the problem 
  
sol = optimize(constr, gamma, opts);

% -- Check the constraints

disp('1) *Constraints check*')
if all(check(constr) > -fake_zero)
disp(['    ' 'Constraints are satisfied'])
else
    disp('Constraints are not satisfied')
end

% -- Extract the results

W = double(W);
X = double(X);
K = X / W;
gamma = double(gamma);

% -- Check stability and convergence rate 

real_eig = real(eig(Anew+B*K)); %Convergence rate --> real part of minimum eigenvalue

disp('2) *Stability test*')
disp(['    ' '*Closed-loop eigenvalues real parts: '])
disp(['        ', mat2str(real_eig')]);
if all(real_eig < 0)
    disp(['    ' 'all eigenvalues have negative real part --> system is stable'])
else
    disp(['    ' 'at least one eigenvalue has positive real part --> system is unstable'])
end

conv_rate = real(max(real_eig));
disp(['Convergence rate: ', num2str(conv_rate)]);

disp('3) *Gain matrix K:')
disp(['    ', mat2str(K)]);
disp('4) *Minimized L2-gain:')
disp(['    ', mat2str(gamma)]);

%% SIMULATION AND PLOTTING

t_f = 4;
max_step = 1e-3;
x0 = [0; 0; 0; 0];
Amplitude = 0.05;
w_sine = 80;
phi_sine = 0;
w_cosine = 30;
phi_cosine = pi/3;

sim("simproject.slx");
t_sim = w.Time; 
w = w.Data;
u = u.Data;

figure
plot(t_sim, w, 'LineWidth',2,'Color',[0.9290 0.6940 0.1250]) 
title('Evolution of the disturbance signal w')
xlabel('t (s)') 
ylabel('Disturbance signal w') 
figure
hold on
plot(t_sim, z_OL, 'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
plot(t_sim, z, 'LineWidth',2,'Color',[0 0.4470 0.7410])
title('Evolution of the performance output z')
xlabel('t (s)') 
ylabel('Performace output z') 
legend('Open-loop case','Closed-loop case','Location','best')
hold off
figure
plot(t_sim, u, 'LineWidth',2,'Color',[0.9290 0.6940 0.1250]) 
title('Closed-loop case: evolution of the control input u')
xlabel('t (s)') 
ylabel('Control input u') 