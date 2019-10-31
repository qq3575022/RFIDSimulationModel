clc, clear, close all

global sim
sim = 1; % sim = 1 for 1D; sim = 2 for 2D

if sim == 1
    clc, clear, close all
    
    load('PropGp1DPositionData.mat')
    load('PropGp1DVelocityAccelerationData.mat')
    
    ta = 0.41;
    tb = 1.41;

    pp = Position_1D(Time >= ta & Time <= tb)'; %-0.99
    vv = MeasVel1D(Time >= ta & Time <= tb)';
    aa = MeasAccel1D(Time >= ta & Time <= tb)';

    td = Time(Time >= ta & Time <= tb)'-ta;
    T = mean(diff(td));
    
    [PP, VV, AA, AStates] = groundtruth1D(td);
    
    % Parameters of tag
    Gt = 14.62;     % tag's antenna gain
    X = 0.5;        % polarization mismatch
    M = 0.25;       % load modulation factor of the tag
    f = 5.8*10^9;
    
    % Parameters of reader
    PT = 1;         % reader's transmitted power
    GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
    GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
    R = 15;

    % Channel noise standard deviation
    sigma = 0.0001; 

    % position of reader
    x0 = [0, 0];
    
    h = 1e-6; tc = 0:h:1-h;
    
    % phase concatenation and amplitude filter
    global l;  l = 0;    k = 1;     
    
    % set state vector
    z_prev = NaN(2,length(td)-1);          % x y position at time step t-1
    z = NaN(2,length(td)-1);               % x y position at time step t

    z_prev(1,:) = PP(1,1:end-1);         z(1,:) = PP(1,2:end);            % x coordinate
    z_prev(2,:) = zeros(1,length(td)-1); z(2,:) = zeros(1,length(td)-1);  % y coordinate
    
    % get simulated results of RSS and phase 
    H         = NaN(1,length(td));         % filtered RSS
    phi_conca = NaN(1,length(td));         % concatenated phase
    phi_mod   = NaN(1,length(td));         % modular phase
    
    % get simulated results of solved r and rdot
    r_sim    = NaN(1,length(td));          % r distance between tag and reader
    rdot_sim = NaN(1,length(td));          % r dot between tag and reader
    
    % get filter buffer
    diff = NaN(1,length(td)-1);
    
    % Get simulation and channel noise filtering
    for k = 1:1:length(td)-1  
        [H(k+1), phi_conca(k+1),phi_mod(k+1),r_sim(k+1),rdot_sim(k+1), diff(k)] = noisysim(x0,f,Gt,M,X,PT,GT,GR,R,sigma,0,k,z,z_prev,T);
    end
    
    H(1) = H(2); phi_conca(1) = phi_conca(2); phi_mod(1) = phi_mod(2); r_sim(1) = r_sim(2); rdot_sim(1) = rdot_sim(2);
 
    errorsim = r_sim-PP;
    error_sim = [mean(errorsim), var(errorsim),sqrt(var(errorsim))]
    cali_errorsim = errorsim;
    
    errormeas = pp-PP;
    error_meas = [mean(errormeas), var(errormeas),sqrt(var(errormeas))]
    cali_errormeas = errormeas;
    
    figure
    subplot(2,1,1),plot(td,H,'LineWidth',2);       title('Simulated H in 1D Motion');    ylabel('Magnitude [V]');xlabel('t [s]')
    subplot(2,1,2),plot(td,phi_mod,'LineWidth',2); title('Simulated phase in 1D Motion');ylabel('phi [rad]');xlabel('t [s]')

    figure
    subplot(2,1,1),
    plot(td,pp,'-',td,r_sim,'LineWidth',2); hold on; plot(td,PP,'LineWidth',1);ylabel('Distance [m]');title('1D Position');
    legend('Measurement','Simulation','Ground truth')

    subplot(2,1,2),
    plot(td,vv,'-',td,rdot_sim,'LineWidth',2); hold on; plot(td,VV,'LineWidth',1);ylabel('Velocity [m/s]');title('1D Velocity');
    legend('Measurement','Simulation','Ground truth')
    xlabel('t [s]')

    figure
    subplot(2,1,1),plot(td,abs(cali_errorsim),'LineWidth',2);xlabel('t [s]');ylabel('|e| [m]');title('Absolute Error of Simulated Position in 1D Motion')
    subplot(2,1,2),plot(td,abs(cali_errormeas),'LineWidth',2);xlabel('t [s]');ylabel('|e| [m]');title('Absolute Error of Measurement Position in 1D Motion')
    

    % --------------------------- For Debug ---------------------------    
%     figure;
%     subplot(3,1,1), plot(td, pp, 'LineWidth',2);hold on; plot(td, PP, 'LineWidth',2);legend('measurement','ground truth'); title('position')
%     subplot(3,1,2), plot(td, vv, 'LineWidth',2);hold on; plot(td, VV, 'LineWidth',2);legend('measurement','ground truth'); title('velocity')
%     subplot(3,1,3), plot(td, aa, 'LineWidth',2);hold on; plot(td, AA, 'LineWidth',2);legend('measurement','ground truth'); title('acceleration')
%   

%     figure
%     subplot(2,1,1),plot(td,phi_mod); title('phi mod');xlabel('t [s]')
%     subplot(2,1,2),plot(td,phi_conca); title('phi concatenated');xlabel('t [s]')
     
else
    
    clc, clear, close all
    load('yVectorData.mat')
    
    offset = 0;
    for n = 2:1:length(psif1)
        if psif1(n) - psif1(n-1) < -2
            offset = 2*pi;
        end
        psif1(n) = psif1(n) + offset - pi;  
    end
    
    psif1(1) = psif1(1) - pi;

    h = 1e-4;           tc = 0:h:1.9099;
    T = mean(diff(tf)); td = 0:T:1.9099;

    % Predefined coordinates of the moving tag
    [W_A, W_V, W, XX, rr] = groundtruth2D(td);
    
    W_VN = W_V + 0.1*rand(1,length(W_V)); WN = W + 0.1*rand(1,length(W));
    
    ax_noisy = rr(7,:) + 0.2*randn(size(rr(7,:)));ay_noisy = rr(8,:) + 0.2*randn(size(rr(8,:)));
    
    figure
    subplot(5,1,1),plot(td,W_A,'LineWidth',2),title('Angular acceleration $\ddot{\psi}$','interpreter','latex','fontweight','bold'),ylabel('$\ddot{\psi}$ [rad/$s^2$]','interpreter','latex'),xlabel('t [s]'), grid%,xlim([0,td(end)])
    subplot(5,1,2),plot(td,W_V,'LineWidth',2),title('Angular velocity $\dot{\psi}$','interpreter','latex','fontweight','bold'),ylabel('$\dot{\psi}$ [rad/$s$]','interpreter','latex'),xlabel('t [s]'), grid%,xlim([0,td(end)])
    subplot(5,1,3),plot(td,W,'LineWidth',2),  title('Orientation $\psi$','interpreter','latex','fontweight','bold'),     ylabel('$\psi$ [rad/s]','interpreter','latex'),     xlabel('t [s]'), grid%,xlim([0,td(end)])
    subplot(5,1,4),plot(td,XX(1,:),'LineWidth',2), title('Position $x$','interpreter','latex','fontweight','bold'), xlabel('t [s]'),ylabel('x[m]','interpreter','latex'), grid
    subplot(5,1,5),plot(td,XX(4,:),'LineWidth',2), title('Position $y$','interpreter','latex','fontweight','bold'),xlabel('t [s]'),ylabel('y[m]','interpreter','latex'), grid
    
    figure
    plot(td,W_V,'LineWidth',1.8),hold on,plot(td,0.85*wf,'LineWidth',1.8),    title('Angular acceleration $\ddot{\psi}$','interpreter','latex','fontweight','bold'),ylabel('$\ddot{\psi}$ [rad/$s^2$]','interpreter','latex'),xlabel('t [s]'), grid%,xlim([0,td(end)])

    % Parameters of tag
    Gt = 14.62;    % tag's antenna gain
    X = 0.85;      % polarization mismatch
    M = 4;         % load modulation factor of the tag
    f = 5.8*10^9;

    % Parameters of reader
    PT = 1;         % reader's transmitted power
    GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
    GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
    R = 15;
    
    % Channel noise error covariance
    sigma = 0.0002; 

    % position of readers
    x1 = [-0.05, 1.5];x2 = [2, 3];x3 = [2.7, 0.05];
    
    % phase cconcatenation
    global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 

    % Get RSS and phase from each reader observing the moving tag
    z = NaN(2,length(td)-1); z_prev = NaN(2,length(td)-1);

    z_prev(1,:) = XX(1,1:end-1); z(1,:) = XX(1,2:end);% x coordinate
    z_prev(2,:) = XX(4,1:end-1); z(2,:) = XX(4,2:end);% y coordinate

    H1 = NaN(1,length(td));phi1 = NaN(1,length(td)); phi_mu1 = NaN(1,length(td));
    H2 = NaN(1,length(td));phi2 = NaN(1,length(td)); phi_mu2 = NaN(1,length(td));
    H3 = NaN(1,length(td));phi3 = NaN(1,length(td)); phi_mu3 = NaN(1,length(td));
    
    r_sim1 = NaN(1,length(td));rdot_sim1 = NaN(1,length(td));diff1 = NaN(1,length(td));
    r_sim2 = NaN(1,length(td));rdot_sim2 = NaN(1,length(td));diff2 = NaN(1,length(td));
    r_sim3 = NaN(1,length(td));rdot_sim3 = NaN(1,length(td));diff3 = NaN(1,length(td));
  
    for k = 1:1:length(td)-1  
    [H1(k+1),phi1(k+1),phi_mu1(k+1),r_sim1(k+1),rdot_sim1(k+1),diff1(k+1)] = noisysim(x1,f,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,T);
    [H2(k+1),phi2(k+1),phi_mu2(k+1),r_sim2(k+1),rdot_sim2(k+1),diff2(k+1)] = noisysim(x2,f,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,T);
    [H3(k+1),phi3(k+1),phi_mu3(k+1),r_sim3(k+1),rdot_sim3(k+1),diff3(k+1)] = noisysim(x3,f,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,T);   
    end
    
    H1(1) = H1(2);phi_mu1(1) = phi_mu1(2);r_sim1(1) = r_sim1(2);rdot_sim1(1)=0;
    H2(1) = H2(2);phi_mu2(1) = phi_mu2(2);r_sim2(1) = r_sim2(2);rdot_sim2(1)=0;
    H3(1) = H3(2);phi_mu3(1) = phi_mu3(2);r_sim3(1) = r_sim3(2);rdot_sim3(1)=0;
    
    % Calculate x, y from r1 r2 r3
    xy_meas = rtoxy(r_sim1, r_sim2, r_sim3, T);
    
    errormea1 = r1f'-rr(1,:);   errormea2 = r2f'-rr(3,:);   errormea3 = r3f'-rr(5,:);
    errorsim1 = r_sim1-rr(1,:); errorsim2 = r_sim2-rr(3,:); errorsim3 = r_sim3-rr(5,:);
   
    error_matrics_meas1 = [mean(errormea1), var(errormea1),sqrt(var(errormea1))]
    error_matrics_meas2 = [mean(errormea2), var(errormea2),sqrt(var(errormea2))]
    error_matrics_meas3 = [mean(errormea3), var(errormea3),sqrt(var(errormea3))]
    
    error_matrics_sim1 = [mean(errorsim1), var(errorsim1),sqrt(var(errorsim1))]
    error_matrics_sim2 = [mean(errorsim2), var(errorsim2),sqrt(var(errorsim2))]
    error_matrics_sim3 = [mean(errorsim3), var(errorsim3),sqrt(var(errorsim3))]
    
    cali_errormea1 = errormea1; cali_errormea2 = errormea2; cali_errormea3 = errormea3 ;
    cali_errorsim1 = errorsim1; cali_errorsim2 = errorsim2; cali_errorsim3 = errorsim3;

    figure
    subplot(3,1,1),plot(td,H1,'LineWidth',2);title('Simulated H from Reader $\#1$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
    subplot(3,1,2),plot(td,H2,'LineWidth',2);title('Simulated H from Reader $\#2$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
    subplot(3,1,3),plot(td,H3,'LineWidth',2);title('Simulated H from Reader $\#3$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')

    figure
    subplot(3,1,1),plot(td,phi_mu1,'LineWidth',2); title('Simulated Phase from Reader $\#1$ in 2D','interpreter','latex');ylabel('phi [rad]');xlabel('t [s]')
    subplot(3,1,2),plot(td,phi_mu2,'LineWidth',2); title('Simulated Phase from Reader $\#2$ in 2D','interpreter','latex');ylabel('phi [rad]');xlabel('t [s]')
    subplot(3,1,3),plot(td,phi_mu3,'LineWidth',2); title('Simulated Phase from Reader $\#3$ in 2D','interpreter','latex');ylabel('phi [rad]');xlabel('t [s]')

    figure
    subplot(3,1,1),plot(td,r1f,'-',td,r_sim1,'LineWidth',2);hold on;plot(td,rr(1,:),'LineWidth',1);legend('Measurement','Simulated','Ground truth');title('2D Radial Distance $R_1$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
    subplot(3,1,2),plot(td,r2f,'-',td,r_sim2,'LineWidth',2);hold on;plot(td,rr(3,:),'LineWidth',1);legend('Measurement','Simulated','Ground truth');title('2D Radial Distance $R_2$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
    subplot(3,1,3),plot(td,r3f,'-',td,r_sim3,'LineWidth',2);hold on;plot(td,rr(5,:),'LineWidth',1);legend('Measurement','Simulated','Ground truth');title('2D Radial Distance $R_3$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')

    figure
    subplot(3,1,1),plot(td,r1dot_e,'-',td,rdot_sim1,'LineWidth',2);hold on;plot(td,rr(2,:),'LineWidth',1);legend('Measurement','Simulated','Ground truth');title('2D Radial Velocity $\dot R_1$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
    subplot(3,1,2),plot(td,r2dot_e,'-',td,rdot_sim2,'LineWidth',2);hold on;plot(td,rr(4,:),'LineWidth',1);legend('Measurement','Simulated','Ground truth');title('2D Radial Velocity $\dot R_2$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
    subplot(3,1,3),plot(td,r3dot_e,'-',td,rdot_sim3,'LineWidth',2);hold on;plot(td,rr(6,:),'LineWidth',1);legend('Measurement','Simulated','Ground truth');title('2D Radial Velocity $\dot R_3$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')

    figure
    subplot(4,1,1),plot(td,0.8*(ax1+0.6),'LineWidth',2);hold on; plot(td,ax_noisy,'LineWidth',2);hold on; plot(td,rr(7,:),'LineWidth',2);legend('Measurement','Simulated','Ground truth');title('Acceleration along $x^B$ axis','interpreter','latex');ylabel('acceleration [m/s^2]');xlabel('t [s]')
    subplot(4,1,2),plot(td,0.8*(ay1+0.15),'LineWidth',2);hold on; plot(td,ay_noisy,'LineWidth',2);hold on; plot(td,rr(8,:),'LineWidth',2);legend('Measurement','Simulated','Ground truth');title('Acceleration along $y^B$ axis','interpreter','latex');ylabel('acceleration [m/s^2]');xlabel('t [s]')
    subplot(4,1,3),plot(td, 0.85*wf, td, W_VN, td, W_V, 'LineWidth',2),legend('Measurement','Simulated','Ground truth'),title('Angular velocity $\dot{\psi}$','interpreter','latex','fontweight','bold'),ylabel('gyroscope output $\dot{\psi}$ [rad/$s$]','interpreter','latex'),xlabel('t [s]')
    subplot(4,1,4),plot(td, psif1, td, WN, td, W, 'LineWidth',2),legend('Measurement','Simulated','Ground truth'),title('Orientation $\psi$','interpreter','latex','fontweight','bold'),ylabel('magnetometer output $\psi$ [rad]','interpreter','latex'),xlabel('t [s]')
    
    figure
    subplot(6,1,1),plot(td,abs(cali_errormea1),'LineWidth',2);title('Absolute error r1');xlabel('t [s]');ylabel('|e| [m]');legend('meas error r1');title('Measurement Absolute Radial Error from Reader $\#1$','interpreter','latex')
    subplot(6,1,2),plot(td,abs(cali_errorsim1),'LineWidth',2);title('Absolute error r1');xlabel('t [s]');ylabel('|e| [m]');legend('sim error r1');title('Simulated Absolute Radial Error from Reader $\#1$','interpreter','latex')

    subplot(6,1,3),plot(td,abs(cali_errormea2),'LineWidth',2);title('Absolute error r2');xlabel('t [s]');ylabel('|e| [m]');legend('meas error r2');title('Measurement Absolute Radial Error from Reader $\#2$','interpreter','latex')
    subplot(6,1,4),plot(td,abs(cali_errorsim2),'LineWidth',2);title('Absolute error r2');xlabel('t [s]');ylabel('|e| [m]');legend('sim error r2');title('Simulated Absolute Radial Error from Reader $\#2$','interpreter','latex')

    subplot(6,1,5),plot(td,abs(cali_errormea3),'LineWidth',2);title('Absolute error r3');xlabel('t [s]');ylabel('|e| [m]');legend('meas error r3');title('Measurement Absolute Radial Error from Reader $\#3$','interpreter','latex')
    subplot(6,1,6),plot(td,abs(cali_errorsim3),'LineWidth',2);title('Absolute error r3');xlabel('t [s]');ylabel('|e| [m]');legend('sim error r3');title('Simulated Absolute Radial Error from Reader $\#3$','interpreter','latex')


    % Compare measurement position, velocity and accelerqtion with ground truth
%     figure
%     subplot(3,1,1),plot(td,r1f,'LineWidth',2);hold on;plot(td,rr(1,:),'LineWidth',2);legend('measurement', 'ground truth');title('r1')
%     subplot(3,1,2),plot(td,r2f,'LineWidth',2);hold on;plot(td,rr(3,:),'LineWidth',2);legend('measurement', 'ground truth');title('r2')
%     subplot(3,1,3),plot(td,r3f,'LineWidth',2);hold on;plot(td,rr(5,:),'LineWidth',2);legend('measurement', 'ground truth');title('r3')
%     
%     figure
%     subplot(3,1,1),plot(td,r1dot_e,'LineWidth',2);hold on;plot(td,rr(2,:),'LineWidth',2);legend('measurement', 'ground truth');title('r1 dot')
%     subplot(3,1,2),plot(td,r2dot_e,'LineWidth',2);hold on;plot(td,rr(4,:),'LineWidth',2);legend('measurement', 'ground truth');title('r2 dot')
%     subplot(3,1,3),plot(td,r3dot_e,'LineWidth',2);hold on;plot(td,rr(6,:),'LineWidth',2);legend('measurement', 'ground truth');title('r3 dot')
%     
%     subplot(2,1,2),
%     plot(td,rdot_sim1,'.',td,rr(2,:),'LineWidth',2); legend('simulation', 'ground truth');hold on;
%     plot(td,rdot_sim2,'.',td,rr(4,:),'LineWidth',2); legend('simulation', 'ground truth');hold on;
%     plot(td,rdot_sim3,'.',td,rr(6,:),'LineWidth',2); legend('simulation', 'ground truth');title('Velocity');
%     legend('R1 dot','Actual R1 dot','R2 dot','Actual R2 dot','R3 dot','Actual R3 dot')
%     xlabel('t [s]')
% 
%     figure
%     subplot(2,1,1),plot(td,xy_meas(1,:),'LineWidth',2);hold on; plot(td,XX(1,:),'LineWidth',2);legend('meas','groundtruth');title('position x')
%     subplot(2,1,2),plot(td,xy_meas(4,:),'LineWidth',2);hold on; plot(td,XX(4,:),'LineWidth',2);legend('meas','groundtruth');title('position y')
%     
%     figure
%     subplot(2,1,1),plot(td,xy_meas(2,:),'LineWidth',2);hold on; plot(td,XX(2,:),'LineWidth',2);legend('meas','groundtruth');title('velocity x')
%     subplot(2,1,2),plot(td,xy_meas(5,:),'LineWidth',2);hold on; plot(td,XX(5,:),'LineWidth',2);legend('meas','groundtruth');title('velocity y')
%     
%     figure
%     subplot(2,1,1),plot(td,xy_meas(3,:),'LineWidth',2);hold on; plot(td,XX(3,:),'LineWidth',2);legend('meas','groundtruth');title('acceleration x')
%     subplot(2,1,2),plot(td,xy_meas(6,:),'LineWidth',2);hold on; plot(td,XX(6,:),'LineWidth',2);legend('meas','groundtruth');title('acceleration y')
%     
%     figure
%     subplot(3,1,1),plot(td,phi3,'LineWidth',2); title('phi 1 concatenated');
%     subplot(3,1,2),plot(td,phi_mu3,'LineWidth',2); title('phi 1 mod');
%     subplot(3,1,3),plot(td,diff3,'LineWidth',2); title('phi 1 diff');
%     xlabel('t [s]')

end


