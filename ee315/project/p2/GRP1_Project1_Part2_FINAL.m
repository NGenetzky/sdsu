% EE315 Project 1 - Part 2 
% Nathan G and Patrick S
% Revised: 03/15/2015
% Version: 12


clearvars 
close all hidden
clc
format compact
warning ('off','all')
%% Options for Output of Running Script
MotorNum = 12;

C1.txtTitle = 'C1.Open-Loop Step Response – neglecting gravity(MATLAB)';
    C1.plotGraphs = 1;
    C1.showPerformanceInfo = 1; 
C2.txtTitle = 'C2.Open-Loop Step Response – neglecting gravity (SIMULINK)';
    C2.runSimulink = 1;
    C2.plotGraphs = 1;
    C2.showPerformanceInfo = 1;
C3.txtTitle = 'C3.Open-Loop Matlab vs Simulink (Neglecting Gravity)';
    C3.plotGraphs = 1; 
D.txtTitle = 'D.Open-Loop Step Response – including effect of gravity (SIMULINK)';
    D.runSimulink = 1;
    D.plotGraphs = 1;
    D.showPerformanceInfo = 1;
D1.txtTitle = 'D1.Compare Open-Loop Simulink with vs without Gravity';
    D1.plotGraphs = 1; 
E1.txtTitle = 'E1.Closed-Loop Step Response neglecting gravity (MATLAB)';
    E1.plotGraphs = 1;
    E1.showPerformanceInfo = 1; 
E2.txtTitle = 'E2.Closed-Loop Step Response including effect of gravity (SIMULINK)';
    E2.runSimulink = 1;
    E2.plotGraphs = 1;
    E2.showPerformanceInfo = 1;

Figures.Units.distance = '[m]';
Figures.Units.time = '[ms]';
Figures.Units.percent = '%';
Figures.Units.velocity = '[m/s]';
Figures.Units.torque = '[N*m]';
Figures.Units.current = '[A]';

Tol = 0.02;             % tolerance level
RH.Mode = 2;            % RH mode selection (1 = simple symbolic (b0,a3,a2,a1,a0), 2 = expanded symbolic (J,Deq,Ra,La,...)            
RH.showTable = 1;     % show RH table & k inequalites (1 = yes, !1 = no)

clean = 1;
SavePNGImages = 1;
OutputFileDir = ['EE315+Part2+Motor',num2str(MotorNum)];
    if ~exist(OutputFileDir, 'dir')
        mkdir(OutputFileDir);
    end
%% Header
disp('--------------------------------------------------');
disp('EE315 Linear Control Systems Project 1 Part 2')
disp('Written By : Nathan G and Patrick S')
disp('EE315_Project1_Part2 Date: 03/16/2015 ')
disp('Wait until the console reappears before any action');
disp('--------------------------------------------------');
%% Dr. Hietpas Variables
if(MotorNum == 4)
    Jm = 6.7; % Kg*m^2
    Dm = 0.117; % N*m*s*(rad^-1)
    kt = 4.4;% N*m*(A^-1)
    kb = 4.6; % V*(rad^-1)*s
    %kc = inf; % N*(m^-1)
    Ra = 29 * 10^(-3); % Ohms 
    La = 0.61* 10^(-3); % H
    UN = 440; % Vdc
end
if(MotorNum == 11)
    Jm = 6.5; % Kg*m^2
    Dm = 0.117; % N*m*s*(rad^-1)
    kt = 3250/747; % N*m*(A^-1) % Numerical Value (4.350736278447122) From Dr. Hietpas (Need To Discuss) was((4.936+5.5)/2)
    kb = kt; % V*(rad^-1)*s
    %kc = inf; % N*(m^-1)
    Ra = 31.1 * 10^(-3); % Ohms 
    La = 0.51* 10^(-3); % H % Numerical Value From Dr. Hietpas was(0.92 mH) ? Double Checked with datasheet = (0.51 mH)
    UN = 440; % Vdc
end
if(MotorNum == 12)
    Jm = 11;        % Kg*m^2
    Dm = 0.117;     % N*m*s*(rad^-1)
    Ra = 54.8 *10^-3;      % Ohms
    La = 0.92 *10^-3;      % H
    UN = 470;       % Vdc
    kt = 3285/675;      % V*(rad^-1)*s      %kt=T/IN; IN = 675, T=3285
    kb = (UN-Ra*675)/(825*(2*pi/60));% N*(m^-1) %kb=(UN-Ra*IN)/(n(2pi/60)) % n=825
end
    % Dr. Hietpas Variables (For all motors)
    Jg = 250;       % Kg*m^2
    Me = 7500;      % Kg
    r = 1;          % m
    Dg = 75.8;      % N*m*s*(rad^-1)
    De = 151.6;     % N*s*(m^-1) 
    Nm = 60;        % Gear Teeth 
    Ng = 1440;      % Gear Teeth

%Other Variables
    
    Deq = Dm + 4*Dg*(Nm/Ng)^2+2*De*r^2*(Nm/Ng)^2;
    J = Jm + 4*Jg*(Nm/Ng)^2+Me*r^2*(Nm/Ng)^2;
    Ag = 9.82;
    Tg= (Nm/Ng)*r*Ag*Me;
    Height = 200*0.3048; % [m] = 200ft*(0.3048 m/ft)
    % K %Defined Later
    %SimulinkInputSelector % Defined later
    %Tg  %Defined later.

%% Determine Transfer Functions for Gx and Gy
%Transfer function Gx
    Gx.b0 = kt*r*(Nm/Ng)*1/(La*J);
    Gx.a0 = 0;
    Gx.a1 = (Deq*Ra+kt*kb)/(La*J);
    Gx.a2 = (Ra*J+Deq*La)/(La*J);
    Gx.a3 = 1;
    Gx.num = Gx.b0;
    Gx.den = [Gx.a3 Gx.a2 Gx.a1 Gx.a0];
    Gx.TF = tf(Gx.num,Gx.den);
    Gx.TF.InputUnit = 'position';
    Gx.TF.OutputUnit = 'voltage';

%Transfer function Gy
    Gv.b0 = Gx.b0;
    Gv.a0 = Gx.a1;
    Gv.a1 = Gx.a2;
    Gv.a2 = Gx.a3;
    Gv.num = [Gv.b0];
    Gv.den = [Gv.a2 Gv.a1 Gv.a0];
    Gv.TF = tf(Gv.num,Gv.den);
    Gv.TF.InputUnit = 'velocity';
    Gv.TF.OutputUnit = 'voltage';
    
    Gv.sigma_d = Gv.a1 / 2;     %Damping Coefficient [Hz] 
    Gv.tau_d = 1/Gv.sigma_d;    % [sec]
    Gv.omega_n = sqrt(Gv.a0);   % Natural Frequency [Hz]
    

                                       % symbolic variables (required for all modes)
%% RH_real table
syms k e 
delta = [Gx.a3 Gx.a2 Gx.a1 Gx.a0+k*Gx.b0];
epsilon = e;

coeff = length(delta);                          % find highest order of coeffs
RH_real = sym(zeros(coeff,round(coeff/2)));     % initialize array 

for i=1:coeff,
	RH_real(2-mod(i,2),round(i/2))=delta(i);    % create 1st and 2nd rows of a coeffs
end

rows=coeff-2;                                   % # of rows to by solved (b,c,d,...)
index=zeros(rows,1);                            % create a column of zeros per row 

for i=1:rows,
	index(rows-i+1)=round(i/2);                 % create index vector from bottom to top
end

for i=3:coeff,                                  % go from 3rd row to last (b-z)
	if(all(RH_real(i-1,:)==0)),                 % special case for row of zeros
			a=coeff-i+2;                        % order of equation
			b=round(a/2)-mod(a,2)+1;            % number of coefficients
            temp1=RH_real(i-2,1:b);             % find polynomial
			temp2=a:-2:0;                       % create the orders for the polynomial
			RH_real(i-1,1:b)=temp1.*temp2;      % derivative of auxiliary
	elseif(RH_real(i-1,1)==0),                  % special case if first element in row is zero
			RH_real(i-1,1)=epsilon;             % speacial case if divide by 0 (replace 0 with epsilon)
    end
                                                % compute the Routh array elements
	for j=1:index(i-2),	
		RH_real(i,j)=-det([RH_real(i-2,1) RH_real(i-2,j+1);RH_real(i-1,1) RH_real(i-1,j+1)])/RH_real(i-1,1);
    end
end
 
% analyze results
RH.k_real = sym(zeros(length(delta),1));         % create array of zeros for inequalities

for i = 2:length(delta)
        RH.k_real(i) = solve(RH_real(i,1),k);    % solve for k (creates inequalities)
end

Kmax = double(RH.k_real(3));                     % determine max value of K (in RH.Mode = 3)
K = Kmax/2;
%% RH_sym table
% save values for reset
RH.Nm_temp = Nm;
RH.Ng_temp = Ng;
RH.r_temp = r;
RH.La_temp = La;
RH.Ra_temp = Ra;
RH.J_temp = J;
RH.Deq_temp = Deq;
RH.kb_temp = kb;
RH.kt_temp = kt;

    if      RH.Mode == 1;                   % simplified symbolic form
        syms a3 a2 a1 a0 b0; 
    elseif  RH.Mode == 2;                   % expanded symbolic form
        syms a3 a2 a1 a0 b0; 
        syms Nm Ng r La Ra J Deq kb kt;  
        b0 = kt*r*(Nm/Ng)*1/(La*J);
        a0 = 0;
        a1 = (Deq*Ra+kt*kb)/(La*J);
        a2 = (Ra*J+Deq*La)/(La*J);
        a3 = 1;
    end
% setup 
delta = [a3 a2 a1 a0+k*b0];
epsilon = e;

coeff = length(delta);                   % find highest order of coeffs
RH_sym = sym(zeros(coeff,round(coeff/2)));    % initialize array 

for i=1:coeff,
	RH_sym(2-mod(i,2),round(i/2))=delta(i);  % create 1st and 2nd rows of a coeffs
end

rows=coeff-2;                           % # of rows to by solved (b,c,d,...)
index=zeros(rows,1);                    % create a column of zeros per row 

for i=1:rows,
	index(rows-i+1)=round(i/2);         % create index vector from bottom to top
end

for i=3:coeff,                          % go from 3rd row to last (b-z)
	if(all(RH_sym(i-1,:)==0)),              % special case for row of zeros
			a=coeff-i+2;                % order of equation
			b=round(a/2)-mod(a,2)+1;    % number of coefficients
            temp1=RH_sym(i-2,1:b);          % find polynomial
			temp2=a:-2:0;               % create the orders for the polynomial
			RH_sym(i-1,1:b)=temp1.*temp2;	% derivative of auxiliary
	elseif(RH_sym(i-1,1)==0),               % special case if first element in row is zero
			RH_sym(i-1,1)=epsilon;          % speacial case if divide by 0 (replace 0 with epsilon)
    end
                                        % compute the Routh array elements
	for j=1:index(i-2),	
		RH_sym(i,j)=-det([RH_sym(i-2,1) RH_sym(i-2,j+1);RH_sym(i-1,1) RH_sym(i-1,j+1)])/RH_sym(i-1,1);
    end
end

% analyze results
RH.k_sym = sym(zeros(length(delta),1));     % create array of zeros for inequalities

for i = 2:length(delta)
        RH.k_sym(i) = solve(RH_sym(i,1),k);     % solve for k (creates inequalities)
end

 %display results 
if RH.showTable == 1;
    RH_Table    = RH_sym
    k_equ       = k > RH.k_sym(2:length(delta),1)
    Kmax
    K
end

clearvars Nm Ng r La Ra J Deq kb kt

% reset variables to values
Nm = RH.Nm_temp;
Ng = RH.Ng_temp;
r = RH.r_temp;
La = RH.La_temp;
Ra = RH.Ra_temp;
J = RH.J_temp;
Deq = RH.Deq_temp;
kb = RH.kb_temp;
kt = RH.kt_temp;
%% Closed Loop Transfer Functions
E1.T.K = K;
E1.T.Gx = Gx.TF;
E1.T.FWD = K*Gx.TF;
E1.T.H = tf([1],[1]);
E1.T.TF = feedback(E1.T.FWD,E1.T.H); 
E1.T.TF.InputUnit = 'position';
E1.T.TF.OutputUnit = 'position';
%% Prepare for simulation and calculation
%Open loop
C1.tau = Gv.tau_d;
C1.Time.Start = 0;
C1.Time.Step = C1.tau/100;
C1.Time.Stop = C1.tau*10-C1.Time.Step;
C1.Time.Vector = [C1.Time.Start : C1.Time.Step : C1.Time.Stop];
%Closed loop
[E1.rlocus.R,E1.rlocus.K] = rlocus(E1.T.Gx, K);
E1.tau = abs(real(E1.rlocus.R(2))).^-1; %0.85/10;%Gv.tau_d*10;
E1.Time.Start = 0;
E1.Time.Step = E1.tau/100;
E1.Time.Stop = E1.tau*10-E1.Time.Step;
E1.Time.Vector = [E1.Time.Start : E1.Time.Step : E1.Time.Stop];

C2.tau = C1.tau;
C2.Time.Start = C1.Time.Start;
C2.Time.Step = C1.Time.Step;
C2.Time.Stop = C1.Time.Stop;
C2.Time.Vector = C1.Time.Vector;

C3.tau = C1.tau;
C3.Time.Start = C1.Time.Start;
C3.Time.Step = C1.Time.Step;
C3.Time.Stop = C1.Time.Stop;
C3.Time.Vector = C1.Time.Vector;

D.tau = C1.tau;
D.Time.Start = C1.Time.Start;
D.Time.Step = C1.Time.Step;
D.Time.Stop = C1.Time.Stop;
D.Time.Vector = C1.Time.Vector;

D1.tau = C1.tau;
D1.Time.Start = C1.Time.Start;
D1.Time.Step = C1.Time.Step;
D1.Time.Stop = C1.Time.Stop;
D1.Time.Vector = C1.Time.Vector;

E2.tau = E1.tau;
E2.Time.Start = E1.Time.Start;
E2.Time.Step = E1.Time.Step;
E2.Time.Stop = E1.Time.Stop;
E2.Time.Vector = E1.Time.Vector;


%C1.txtTitle = 'Open-Loop Step Response – neglecting gravity(MATLAB)';
    C1.Gv = Gv;
    C1.figureList=[];
    C1.StepAmplitude = UN;
    C1.ShortTitle = 'C1';
%C2.txtTitle = 'Open-Loop Step Response – neglecting gravity (SIMULINK)';
    C2.Gv = Gv;
    C2.Tg = 0;
    C2.figureList=[];
    C2.SimulinkInputSelector = 1; %Use the Unit Step as input to Gv
    C2.StepAmplitude = UN;
    C2.ShortTitle = 'C2';
%C3.txtTitle = 'C3.Open-Loop Matlab vs Simulink (Neglecting Gravity)';
    C3.ShortTitle = 'C3';
    C3.figureList=[];    
%D.txtTitle =    'Open-Loop Step Response – including effect of gravity (SIMULINK)';
    D.Gv = Gv;
    D.Tg = Tg;
    D.figureList=[];
    D.SimulinkInputSelector = 1; %Use the Unit Step as input to Gv
    D.StepAmplitude = UN;
    D.ShortTitle = 'D';
%D1.txtTitle = 'D1.Compare Open-Loop Simulink with vs without Gravity';
    D1.ShortTitle = 'D1';
    D1.figureList=[];
%E1.txtTitle = 'Closed-Loop Step Response neglecting gravity (MATLAB)';
    E1.figureList=[];
    E1.StepAmplitude = Height;
    E1.ShortTitle = 'E1';
%E2.txtTitle = 'Closed-Loop Step Response including effect of gravity (SIMULINK)';
    E2.Tg =Tg;  
    E2.figureList=[];
    E2.SimulinkInputSelector = 2; %Use the Position as input to Gx
    E2.StepAmplitude = Height;
    E2.ShortTitle = 'E2';
    
%%  C1. Step Response of Gv (using lsim)
    C1.dynamicSim.Input = C1.StepAmplitude*ones(numel(C1.Time.Vector),1);
    C1.dynamicSim.Output = lsim(Gv.TF, C1.dynamicSim.Input,C1.Time.Vector);
    
%% C1. Step Response of Gv (using step)
    C1.Step.Output = step(C1.Gv.TF,C1.Time.Vector,stepDataOptions('StepAmplitude',C1.StepAmplitude));
%Calculate Output from Step Response of Gv (using step)
    %1. % overshoot
    %2. peak speed (Peak Amplitude)
    %3. time to peak 
    %4. steady state speed (Average of Settling Min/Max)
    %5. settling time
    C1.Step.Info = stepinfo(C1.Step.Output, C1.Time.Vector);
    C1.Overshoot = C1.Step.Info.Overshoot;
    C1.PeakVelocity = C1.Step.Info.Peak;
    C1.TimeToPeak = C1.Step.Info.PeakTime;
    C1.SteadyVelocity = C1.Step.Output(numel(C1.Step.Output));
    C1.SettlingTime = C1.Step.Info.SettlingTime;    

%%  E1. Step Response of Closed Loop (using lsim)
    E1.dynamicSim.Input = E1.StepAmplitude*ones(numel(E1.Time.Vector),1);
    E1.dynamicSim.Output = lsim(E1.T.TF, E1.dynamicSim.Input,E1.Time.Vector);
    
%% E1. Step Response of Closed Loop (using step)
    E1.Step.Output = step(E1.T.TF,E1.Time.Vector,stepDataOptions('StepAmplitude',E1.StepAmplitude));
%Calculate Output from Step Response of Gv (using step)
    %1. % overshoot
    %2. peak speed (Peak Amplitude)
    %3. time to peak 
    %4. steady state speed (Average of Settling Min/Max)
    %5. settling time
    E1.Step.Info = stepinfo(E1.Step.Output, E1.Time.Vector);
    E1.Overshoot = E1.Step.Info.Overshoot;
    E1.PeakPosition = E1.Step.Info.Peak;
    E1.TimeToPeak = E1.Step.Info.PeakTime;
    E1.SteadyPosition = E1.Step.Output(numel(E1.Step.Output));
    E1.SettlingTime = E1.Step.Info.SettlingTime;   

%% Create LTI View Plot
%ltiview('step',Gv.TF);

%% Run Simulink
Simulink.ModelName='GRP1_Project1_Part2_Model';
    load_system(Simulink.ModelName);
    Simulink.BlockPaths = find_system(Simulink.ModelName,'Type','Block');
    
if(C2.runSimulink)
    tau = C2.tau;
    Time.Start = C2.Time.Start;
    Time.Step = C2.Time.Step;
    Time.Stop = C2.Time.Stop;
    Time.Vector = C2.Time.Vector;

    Tg = C2.Tg;
    SimulinkInputSelector = C2.SimulinkInputSelector;
    C2.SimOut = sim(Simulink.ModelName, 'ReturnWorkspaceOutputs', 'on');
    % Extract Variables from Simulink
    C2.tout = C2.SimOut.get('tout');
    C2.Va = C2.SimOut.get('Va');
    C2.Ia = C2.SimOut.get('Ia');
    C2.Tm = C2.SimOut.get('Tm');
    C2.am = C2.SimOut.get('am');
    C2.Wm = C2.SimOut.get('Wm');
    C2.Wg = C2.SimOut.get('Wg');
    C2.Vx = C2.SimOut.get('Vx');
    C2.Xx = C2.SimOut.get('Xx');
    C2.Vb = C2.SimOut.get('Vb');
end
if(D.runSimulink)
    tau = D.tau;
    Time.Start = D.Time.Start;
    Time.Step = D.Time.Step;
    Time.Stop = D.Time.Stop;
    Time.Vector = D.Time.Vector;
    
    Tg = D.Tg;
    SimulinkInputSelector = D.SimulinkInputSelector;
    D.SimOut = sim(Simulink.ModelName, 'ReturnWorkspaceOutputs', 'on');
    % Extract Variables from Simulink
    D.tout = D.SimOut.get('tout');
    D.Va = D.SimOut.get('Va');
    D.Ia = D.SimOut.get('Ia');
    D.Tm = D.SimOut.get('Tm');
    D.am = D.SimOut.get('am');
    D.Wm = D.SimOut.get('Wm');
    D.Wg = D.SimOut.get('Wg');
    D.Vx = D.SimOut.get('Vx');
    D.Xx = D.SimOut.get('Xx');
    D.Vb = D.SimOut.get('Vb');
end
if(E2.runSimulink)
    tau = E2.tau;
    Time.Start = E2.Time.Start;
    Time.Step = E2.Time.Step;
    Time.Stop = E2.Time.Stop;
    Time.Vector = E2.Time.Vector;
    
    Tg = E2.Tg;
    SimulinkInputSelector = E2.SimulinkInputSelector;
    E2.SimOut = sim(Simulink.ModelName, 'ReturnWorkspaceOutputs', 'on');
    % Extract Variables from Simulink
    E2.tout = E2.SimOut.get('tout');
    E2.Va = E2.SimOut.get('Va');
    E2.Ia = E2.SimOut.get('Ia');
    E2.Tm = E2.SimOut.get('Tm');
    E2.am = E2.SimOut.get('am');
    E2.Wm = E2.SimOut.get('Wm');
    E2.Wg = E2.SimOut.get('Wg');
    E2.Vx = E2.SimOut.get('Vx');
    E2.Xx = E2.SimOut.get('Xx');
    E2.Vb = E2.SimOut.get('Vb');
end
%% Output Helpers

f_figureTitle = @(X,str) title({X.txtTitle,str});
f_evalAtTime = @(X,t,t0) X(find(t==t0));    % X and T are vectors. Will return X0 at specific time (T0)
f_savePng = @(H,name) print(H,[OutputFileDir,'\',name],'-dpng','-r150');
    %saveas(H,[OutputFilePrefix,name],'png');
f_newFigure = @() figure('units','normalized','outerposition',[0 0 1 1]);

Figures.timeLabel = 't [ms]';

Figures.Va.Title ='Voltage Supplied to Armature';
Figures.Va.Label ='V_a [Volts]'; 
Figures.Ia.Title ='Current through the armature';
Figures.Ia.Label ='I_a [Amp]';
Figures.Tm.Title ='Torque of the motor';
Figures.Tm.Label ='T_m [Nm]';
Figures.am.Title ='Radial Acceleration of the Motor';
Figures.am.Label ='\alpha_m [m/s^2]';
Figures.Wg.Title ='Radial Velocity of the Gears';
Figures.Wg.Label ='\omega_g';
Figures.Vx.Title ='Velocity of Elevator';
Figures.Vx.Label ='V_x [m/s]';
Figures.Xx.Title ='Position of the Elevator';
Figures.Xx.Label ='X_x [m]';
Figures.Vb.Title ='Back EMF from PMDC motor';
Figures.Vb.Label ='V_b [Volts]';
Figures.rlocus.Title = 'Root Locus';

%scrsz = get(groot,'ScreenSize');
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

%% Output for C1 (MATLAB)
    C1.txtOvershoot =       ['1. % overshoot:',num2str(C1.Overshoot),' ',Figures.Units.percent];
    C1.txtPeakVelocity =    ['2. peak speed (Peak Amplitude):',num2str(C1.PeakVelocity),' ',Figures.Units.velocity];
    C1.txtTimeToPeak =      ['3. time to peak:',num2str(C1.TimeToPeak),' ',Figures.Units.time];
    C1.txtSteadyVelocity =  ['4. steady state speed:',num2str(C1.SteadyVelocity),' ',Figures.Units.velocity];
    C1.txtSettlingTime =    ['5. settling time:',num2str(C1.SettlingTime),' ',Figures.Units.time];
if(C1.showPerformanceInfo)
    disp(C1.txtTitle); 
    disp('Performance values for Gv:');
    disp(C1.txtOvershoot); 
    disp(C1.txtPeakVelocity); 
    disp(C1.txtTimeToPeak); 
    disp(C1.txtSteadyVelocity);
    disp(C1.txtSettlingTime);
end
if(C1.plotGraphs)
    C1.figureList = [C1.figureList f_newFigure()];
    %annotations--------------------------------------------------------------------------------------------------------------------
        C1.loc.V_x.PeakOffset = C1.TimeToPeak-C1.TimeToPeak*0.05;
        C1.loc.V_x.Peak = C1.TimeToPeak;
        C1.loc.V_x.Settling = C1.SettlingTime;
        C1.loc.V_y.PeakOffset = C1.PeakVelocity-C1.PeakVelocity*0.03;
        C1.loc.V_y.Peak = C1.PeakVelocity;
        C1.loc.V_y.Steady = C1.SteadyVelocity;
        C1.loc.V_y.Overshoot = (C1.SteadyVelocity+C1.PeakVelocity)/2;
        C1.loc.V_y.UpperTol = C1.SteadyVelocity+C1.SteadyVelocity*Tol;
        C1.loc.V_y.LowerTol = C1.SteadyVelocity-C1.SteadyVelocity*Tol;
        hold on 
        line(C1.loc.V_x.Peak,C1.loc.V_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                    % point for PeakVelocity
        line([C1.loc.V_x.PeakOffset,C1.loc.V_x.Peak],[C1.loc.V_y.PeakOffset,C1.loc.V_y.Peak],'Color',[0 204 0]/255);        % line for PeakVelocity
        line([0,C1.loc.V_x.Peak],[C1.loc.V_y.Peak,C1.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % horizontal line for PeakVelocity
        line([C1.loc.V_x.Peak,C1.loc.V_x.Peak],[0,C1.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % vertical line for PeakVelocity
        line([0,C1.Time.Stop],[C1.loc.V_y.Steady,C1.loc.V_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');               % horizontal line for SteadyVelocity
        line([0,C1.Time.Stop],[C1.loc.V_y.UpperTol,C1.loc.V_y.UpperTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyVelocity UpperTol
        line([0,C1.Time.Stop],[C1.loc.V_y.LowerTol,C1.loc.V_y.LowerTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyVelocity LowerTol
        line([C1.loc.V_x.Settling,C1.loc.V_x.Settling],[0,C1.Step.Output(max(find(abs(((C1.Step.Output-C1.SteadyVelocity)/C1.SteadyVelocity))>Tol))+1)],'Color',[255 153 51]/255,'LineStyle','--');% vertical line for SettlingTime
        line([0.05,0.05],[C1.loc.V_y.Steady,C1.loc.V_y.Peak],'Color',[0 0 0]/255,'LineStyle',':');                      % vertical line for Overshoot
        %------------------------------------------------------------------------------------------------------------------------------- 
        text(C1.loc.V_x.PeakOffset,C1.loc.V_y.PeakOffset,['Peak Velocity = ',num2str(C1.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','top','HorizontalAlignment','right','color',[0 204 0]/255,'FontSize',8); % label for PeakVelocity
        text(0.001,C1.loc.V_y.Peak,['Peak Velocity = ',num2str(C1.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                   % label for horizontal line for PeakVelocity
        text(C1.loc.V_x.Peak,0.1,['  Time to Peak = ',num2str(C1.TimeToPeak),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                      % label for vertical line for PeakVelocity
        text(0.001,C1.loc.V_y.Steady,['Steady State Velocity = ',num2str(C1.SteadyVelocity),' ',Figures.Units.velocity],'VerticalAlignment','top','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);          % label for horizontal line for SteadyVelocity
        text(0.001,C1.loc.V_y.LowerTol,['Steady State Velocity Range: \pm ',num2str(Tol*100),' ',Figures.Units.percent],'VerticalAlignment','top','HorizontalAlignment','left','color',[204 204 0]/255,'FontSize',8);                        % label for horizontal lines for Tolerance Range
        text(C1.loc.V_x.Settling,0.1,['  Settling Time = ',num2str(C1.SettlingTime),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 153 51]/255,'FontSize',8);           % label for vertical line for SettlingTime
        text(0.05,C1.loc.V_y.Overshoot,['  % Overshoot = ',num2str(C1.Overshoot),' ',Figures.Units.percent],'VerticalAlignment','cap','HorizontalAlignment','left','color',[0 0 0]/255,'FontSize',8);                       % label for vertical line for Overshoot
    %-------------------------------------------------------------------------------------------------------------------------------
    plot(C1.Time.Vector,C1.Step.Output,'b');
    f_figureTitle(C1, Figures.Vx.Title);
    xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);  
    hold off
end
%% Calculate Output for C2 (SIMULINK)
if(C2.runSimulink)  
    %Supports negative value by preserving the sign.
    C2.PeakVelocity = max(abs(C2.Vx))*sign(max(abs(C2.Vx)));
    %Finds 'tout' of the index of the values of 'Vx' that are PeakVelocity
    C2.TimeToPeak = C2.tout(find(C2.Vx==C2.PeakVelocity));
    C2.SteadyVelocity = C2.Vx(numel(C2.Vx)-1);
    %Finds 'tout' of the index of the values of 'Vx' that are within 2% (0.02) of Steady State
    C2.SettlingTime = C2.tout(max(find(abs(((C2.Vx-C2.SteadyVelocity)/C2.SteadyVelocity))>Tol))+1);
    C2.Overshoot = (C2.PeakVelocity-C2.SteadyVelocity)/C2.SteadyVelocity *100;

    C2.PeakMotorCurrent = max(C2.Ia);
    C2.TimeToPeakMotorCurrent = C2.tout(find(C2.Ia==C2.PeakMotorCurrent));
    C2.SteadyMotorCurrent = C2.Ia(numel(C2.Ia)-1); % Get last value recorded
    C2.TimeToSteadyMotorCurrent = C2.tout(find(C2.Ia==C2.SteadyMotorCurrent));
    
    C2.PeakMotorTorque = max(C2.Tm);
    C2.TimeToPeakMotorTorque = C2.tout(find(C2.Tm==C2.PeakMotorTorque));
    C2.SteadyMotorTorque = C2.Tm(numel(C2.Tm)-1); % Get last value recorded
    C2.TimeToSteadyMotorTorque = C2.tout(find(C2.Tm==C2.SteadyMotorTorque));

% Prepare Output
    C2.txtOvershoot =       ['1. % overshoot: ',num2str(C2.Overshoot),' ',Figures.Units.percent];
    C2.txtPeakVelocity =    ['2. peak speed (Peak Amplitude): ',num2str(C2.PeakVelocity),' ',Figures.Units.velocity];
    C2.txtTimeToPeak =      ['3. time to peak: ',num2str(C2.TimeToPeak),' ',Figures.Units.time];
    C2.txtSteadyVelocity =  ['4. steady state speed: ',num2str(C2.SteadyVelocity),' ',Figures.Units.velocity];
    C2.txtSettlingTime =    ['5. settling time: ',num2str(C2.SettlingTime),' ',Figures.Units.time];

    C2.txtPeakMotorCurrent =    ['6. peak motor input current: ',num2str(C2.PeakMotorCurrent),' ',Figures.Units.current];
    C2.txtSteadyMotorCurrent =  ['7. steady-state motor input current: ',num2str(C2.SteadyMotorCurrent),' ',Figures.Units.current];
    C2.txtPeakMotorTorque =     ['8. peak motor output torque: ',num2str(C2.PeakMotorTorque),' ',Figures.Units.torque];
    C2.txtSteadyMotorTorque =   ['9. steady-state motor output torque: ',num2str(C2.SteadyMotorTorque),' ',Figures.Units.torque];
end

%% Output for C2 (SIMULINK)
if(C2.runSimulink)
    if(C2.showPerformanceInfo)
        disp(['<<<   ',C2.txtTitle,'   >>>']);
        disp('Performance values:');
        disp(C2.txtOvershoot); 
        disp(C2.txtPeakVelocity); 
        disp(C2.txtTimeToPeak); 
        disp(C2.txtSteadyVelocity);
        disp(C2.txtSettlingTime);

        disp(C2.txtPeakMotorCurrent);
        disp(C2.txtSteadyMotorCurrent);
        disp(C2.txtPeakMotorTorque);
        disp(C2.txtSteadyMotorTorque);
    end
    if(C2.plotGraphs && C2.runSimulink)
        C2.figureList = [C2.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            C2.loc.Tm_x.PeakOffset = C2.TimeToPeakMotorTorque+C2.TimeToPeakMotorTorque*1;
            C2.loc.Tm_x.Peak = C2.TimeToPeakMotorTorque;
            C2.loc.Tm_x.Settling = C2.TimeToSteadyMotorTorque;
            C2.loc.Tm_y.PeakOffset = C2.PeakMotorTorque-C2.PeakMotorTorque*0.03;
            C2.loc.Tm_y.Peak = C2.PeakMotorTorque;
            C2.loc.Tm_y.Steady = C2.SteadyMotorTorque;
            C2.loc.Tm_y.Overshoot = (C2.SteadyMotorTorque+C2.PeakMotorTorque)/2;
            C2.loc.Tm_y.UpperTol = C2.SteadyMotorTorque+C2.SteadyMotorTorque*Tol;
            C2.loc.Tm_y.LowerTol = C2.SteadyMotorTorque-C2.SteadyMotorTorque*Tol;
            hold on 
            line(C2.loc.Tm_x.Peak,C2.loc.Tm_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                          % point for PeakTorque
            line([C2.loc.Tm_x.PeakOffset,C2.loc.Tm_x.Peak],[C2.loc.Tm_y.PeakOffset,C2.loc.Tm_y.Peak],'Color',[0 204 0]/255);        % line for PeakTorque
            line([0,C2.loc.Tm_x.Peak],[C2.loc.Tm_y.Peak,C2.loc.Tm_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                  % horizontal line for PeakTorque
            line([C2.loc.Tm_x.Peak,C2.loc.Tm_x.Peak],[0,C2.loc.Tm_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                  % vertical line for PeakTorque
            line([0,Time.Stop],[C2.loc.Tm_y.Steady,C2.loc.Tm_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');                     % horizontal line for SteadyTorque
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(C2.loc.Tm_x.PeakOffset,C2.loc.Tm_y.PeakOffset,['Peak Torque = ',num2str(C2.PeakMotorTorque),' ',Figures.Units.torque],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);           % label for PeakTorque
            text(C2.loc.Tm_x.Settling,C2.loc.Tm_y.Steady+250,['Steady State Torque = ',num2str(C2.SteadyMotorTorque),' ',Figures.Units.torque],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);    % label for horizontal line for SteadyTorque
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(C2.Time.Vector,C2.Tm);
        f_figureTitle(C2, Figures.Tm.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Tm.Label);
        hold off
        
        C2.figureList = [C2.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            C2.loc.Ia_x.PeakOffset = C2.TimeToPeakMotorCurrent+C2.TimeToPeakMotorCurrent*1;
            C2.loc.Ia_x.Peak = C2.TimeToPeakMotorCurrent;
            C2.loc.Ia_x.Settling = C2.TimeToSteadyMotorCurrent;
            C2.loc.Ia_y.PeakOffset = C2.PeakMotorCurrent-C2.PeakMotorCurrent*0.03;
            C2.loc.Ia_y.Peak = C2.PeakMotorCurrent;
            C2.loc.Ia_y.Steady = C2.SteadyMotorCurrent;
            C2.loc.Ia_y.Overshoot = (C2.SteadyMotorCurrent+C2.PeakMotorCurrent)/2;
            C2.loc.Ia_y.UpperTol = C2.SteadyMotorCurrent+C2.SteadyMotorCurrent*Tol;
            C2.loc.Ia_y.LowerTol = C2.SteadyMotorCurrent-C2.SteadyMotorCurrent*Tol;
            hold on 
            line(C2.loc.Ia_x.Peak,C2.loc.Ia_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                          % point for PeakCurrent
            line([C2.loc.Ia_x.PeakOffset,C2.loc.Ia_x.Peak],[C2.loc.Ia_y.PeakOffset,C2.loc.Ia_y.Peak],'Color',[0 204 0]/255);        % line for PeakCurrent
            line([0,C2.loc.Ia_x.Peak],[C2.loc.Ia_y.Peak,C2.loc.Ia_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                  % horizontal line for PeakCurrent
            line([C2.loc.Ia_x.Peak,C2.loc.Ia_x.Peak],[0,C2.loc.Ia_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                  % vertical line for PeakCurrent
            line([0,Time.Stop],[C2.loc.Ia_y.Steady,C2.loc.Ia_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');                     % horizontal line for SteadyCurrent
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(C2.loc.Ia_x.PeakOffset,C2.loc.Ia_y.PeakOffset,['Peak Current = ',num2str(C2.PeakMotorCurrent),' ',Figures.Units.current],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);           % label for PeakCurrent
            text(C2.loc.Ia_x.Settling,C2.loc.Ia_y.Steady+100,['Steady State Current = ',num2str(C2.SteadyMotorCurrent),' ',Figures.Units.current],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);    % label for horizontal line for SteadyCurrent
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(C2.Time.Vector,C2.Ia);
        f_figureTitle(C2, Figures.Ia.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Ia.Label);
        hold off
        
        C.figureList = [C2.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            C2.loc.V_x.PeakOffset = C2.TimeToPeak-C2.TimeToPeak*0.05;
            C2.loc.V_x.Peak = C2.TimeToPeak;
            C2.loc.V_x.Settling = C2.SettlingTime;
            C2.loc.V_y.PeakOffset = C2.PeakVelocity-C2.PeakVelocity*0.03;
            C2.loc.V_y.Peak = C2.PeakVelocity;
            C2.loc.V_y.Steady = C2.SteadyVelocity;
            C2.loc.V_y.Overshoot = (C2.SteadyVelocity+C2.PeakVelocity)/2;
            C2.loc.V_y.UpperTol = C2.SteadyVelocity+C2.SteadyVelocity*Tol;
            C2.loc.V_y.LowerTol = C2.SteadyVelocity-C2.SteadyVelocity*Tol;
            hold on 
            line(C2.loc.V_x.Peak,C2.loc.V_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                        % point for PeakVelocity
            line([C2.loc.V_x.PeakOffset,C2.loc.V_x.Peak],[C2.loc.V_y.PeakOffset,C2.loc.V_y.Peak],'Color',[0 204 0]/255);        % line for PeakVelocity
            line([0,C2.loc.V_x.Peak],[C2.loc.V_y.Peak,C2.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % horizontal line for PeakVelocity
            line([C2.loc.V_x.Peak,C2.loc.V_x.Peak],[0,C2.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % vertical line for PeakVelocity
            line([0,C2.Time.Stop],[C2.loc.V_y.Steady,C2.loc.V_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');                % horizontal line for SteadyVelocity
            line([0,C2.Time.Stop],[C2.loc.V_y.UpperTol,C2.loc.V_y.UpperTol],'Color',[204 204 0]/255,'LineStyle','--');          % horizontal line for SteadyVelocity UpperTol
            line([0,C2.Time.Stop],[C2.loc.V_y.LowerTol,C2.loc.V_y.LowerTol],'Color',[204 204 0]/255,'LineStyle','--');          % horizontal line for SteadyVelocity LowerTol
            line([C2.loc.V_x.Settling,C2.loc.V_x.Settling],[0,C2.Vx(max(find(abs(((C2.Vx-C2.SteadyVelocity)/C2.SteadyVelocity))>Tol))+1)],'Color',[255 153 51]/255,'LineStyle','--');  % vertical line for SettlingTime
            line([0.05,0.05],[C2.loc.V_y.Steady,C2.loc.V_y.Peak],'Color',[0 0 0]/255,'LineStyle',':');                          % vertical line for Overshoot
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(C2.loc.V_x.PeakOffset,C2.loc.V_y.PeakOffset,['Peak Velocity = ',num2str(C2.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','top','HorizontalAlignment','right','color',[0 204 0]/255,'FontSize',8); % label for PeakVelocity
            text(0.001,C2.loc.V_y.Peak,['Peak Velocity = ',num2str(C2.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                   % label for horizontal line for PeakVelocity
            text(C2.loc.V_x.Peak,0.1,['  Time to Peak = ',num2str(C2.TimeToPeak),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                      % label for vertical line for PeakVelocity
            text(0.001,C2.loc.V_y.Steady,['Steady State Velocity = ',num2str(C2.SteadyVelocity),' ',Figures.Units.velocity],'VerticalAlignment','top','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);          % label for horizontal line for SteadyVelocity
            text(0.001,C2.loc.V_y.LowerTol,['Steady State Velocity Range: \pm ',num2str(Tol*100),' ',Figures.Units.percent],'VerticalAlignment','top','HorizontalAlignment','left','color',[204 204 0]/255,'FontSize',8);                        % label for horizontal lines for Tolerance Range
            text(C2.loc.V_x.Settling,0.1,['  Settling Time = ',num2str(C2.SettlingTime),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 153 51]/255,'FontSize',8);           % label for vertical line for SettlingTime
            text(0.05,C2.loc.V_y.Overshoot,['  % Overshoot = ',num2str(C2.Overshoot),' ',Figures.Units.percent],'VerticalAlignment','cap','HorizontalAlignment','left','color',[0 0 0]/255,'FontSize',8);                       % label for vertical line for Overshoot
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(C2.Time.Vector,C2.Vx);
        f_figureTitle(C2, Figures.Vx.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);
        hold off
    end
end

%% Calculate Output for D (SIMULINK)
if(D.runSimulink)
    %Supports negative value by preserving the sign.
    D.PeakVelocity = max(abs(D.Vx))*sign(max(abs(D.Vx)));
    %Finds 'tout' of the index of the values of 'Vx' that are PeakVelocity
    D.TimeToPeak = D.tout(find(D.Vx==D.PeakVelocity));
    D.SteadyVelocity = D.Vx(numel(D.Vx)-1);
    %Finds 'tout' of the index of the values of 'Vx' that are within 2% (0.02) of Steady State
    D.SettlingTime = D.tout(max(find(abs(((D.Vx-D.SteadyVelocity)/D.SteadyVelocity))>Tol))+1);
    D.Overshoot = (D.PeakVelocity-D.SteadyVelocity)/D.SteadyVelocity *100;

    D.PeakMotorCurrent = max(D.Ia);
    D.TimeToPeakMotorCurrent = D.tout(find(D.Ia==D.PeakMotorCurrent));
    D.SteadyMotorCurrent = D.Ia(numel(D.Ia)-1); % Get last value recorded
    D.TimeToSteadyMotorCurrent = D.tout(find(D.Ia==D.SteadyMotorCurrent));
    
    D.PeakMotorTorque = max(D.Tm);
    D.TimeToPeakMotorTorque = D.tout(find(D.Tm==D.PeakMotorTorque));
    D.SteadyMotorTorque = D.Tm(numel(D.Tm)-1); % Get last value recorded
    D.TimeToSteadyMotorTorque = D.tout(find(D.Tm==D.SteadyMotorTorque));

    % Prepare Output
    D.txtOvershoot =       ['1. % overshoot: ',num2str(D.Overshoot),' ',Figures.Units.percent];
    D.txtPeakVelocity =    ['2. peak speed (Peak Amplitude): ',num2str(D.PeakVelocity),' ',Figures.Units.velocity];
    D.txtTimeToPeak =      ['3. time to peak: ',num2str(D.TimeToPeak),' ',Figures.Units.time];
    D.txtSteadyVelocity =  ['4. steady state speed: ',num2str(D.SteadyVelocity),' ',Figures.Units.velocity];
    D.txtSettlingTime =    ['5. settling time: ',num2str(D.SettlingTime),' ',Figures.Units.time];

    D.txtPeakMotorCurrent =    ['6. peak motor input current: ',num2str(D.PeakMotorCurrent),' ',Figures.Units.current];
    D.txtSteadyMotorCurrent =  ['7. steady-state motor input current: ',num2str(D.SteadyMotorCurrent),' ',Figures.Units.current];
    D.txtPeakMotorTorque =     ['8. peak motor output torque: ',num2str(D.PeakMotorTorque),' ',Figures.Units.torque];
    D.txtSteadyMotorTorque =  ['9. steady-state motor output torque: ',num2str(D.SteadyMotorTorque),' ',Figures.Units.torque];
end
%% Output for D (SIMULINK)
if(D.runSimulink)
    if(D.showPerformanceInfo)
        disp(['<<<   ',D.txtTitle,'   >>>']);
        disp('Performance values:');
        disp(D.txtOvershoot); 
        disp(D.txtPeakVelocity); 
        disp(D.txtTimeToPeak); 
        disp(D.txtSteadyVelocity);
        disp(D.txtSettlingTime);

        disp(D.txtPeakMotorCurrent);
        disp(D.txtSteadyMotorCurrent);
        disp(D.txtPeakMotorTorque);
        disp(D.txtSteadyMotorTorque);
    end
    if(D.plotGraphs  && D.runSimulink)
        D.figureList = [D.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            D.loc.Tm_x.PeakOffset = D.TimeToPeakMotorTorque+D.TimeToPeakMotorTorque*1;
            D.loc.Tm_x.Peak = D.TimeToPeakMotorTorque;
            D.loc.Tm_x.Settling = D.TimeToSteadyMotorTorque;
            D.loc.Tm_y.PeakOffset = D.PeakMotorTorque-D.PeakMotorTorque*0.03;
            D.loc.Tm_y.Peak = D.PeakMotorTorque;
            D.loc.Tm_y.Steady = D.SteadyMotorTorque;
            D.loc.Tm_y.Overshoot = (D.SteadyMotorTorque+D.PeakMotorTorque)/2;
            D.loc.Tm_y.UpperTol = D.SteadyMotorTorque+D.SteadyMotorTorque*Tol;
            D.loc.Tm_y.LowerTol = D.SteadyMotorTorque-D.SteadyMotorTorque*Tol;
            hold on 
            line(D.loc.Tm_x.Peak,D.loc.Tm_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                        % point for PeakTorque
            line([D.loc.Tm_x.PeakOffset,D.loc.Tm_x.Peak],[D.loc.Tm_y.PeakOffset,D.loc.Tm_y.Peak],'Color',[0 204 0]/255);        % line for PeakTorque
            line([0,D.loc.Tm_x.Peak],[D.loc.Tm_y.Peak,D.loc.Tm_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % horizontal line for PeakTorque
            line([D.loc.Tm_x.Peak,D.loc.Tm_x.Peak],[0,D.loc.Tm_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % vertical line for PeakTorque
            line([0,Time.Stop],[D.loc.Tm_y.Steady,D.loc.Tm_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');                   % horizontal line for SteadyTorque
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(D.loc.Tm_x.PeakOffset,D.loc.Tm_y.PeakOffset,['Peak Torque = ',num2str(D.PeakMotorTorque),' ',Figures.Units.torque],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);           % label for PeakTorque
            text(D.loc.Tm_x.Settling,D.loc.Tm_y.Steady,['Steady State Torque = ',num2str(D.SteadyMotorTorque),' ',Figures.Units.torque],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);    % label for horizontal line for SteadyTorque
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(D.Time.Vector,D.Tm);
        f_figureTitle(D, Figures.Tm.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Tm.Label);
        hold off
        
        D.figureList = [D.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            D.loc.Ia_x.PeakOffset = D.TimeToPeakMotorCurrent+D.TimeToPeakMotorCurrent*1;
            D.loc.Ia_x.Peak = D.TimeToPeakMotorCurrent;
            D.loc.Ia_x.Settling = D.TimeToSteadyMotorCurrent;
            D.loc.Ia_y.PeakOffset = D.PeakMotorCurrent-D.PeakMotorCurrent*0.03;
            D.loc.Ia_y.Peak = D.PeakMotorCurrent;
            D.loc.Ia_y.Steady = D.SteadyMotorCurrent;
            D.loc.Ia_y.Overshoot = (D.SteadyMotorCurrent+D.PeakMotorCurrent)/2;
            D.loc.Ia_y.UpperTol = D.SteadyMotorCurrent+D.SteadyMotorCurrent*Tol;
            D.loc.Ia_y.LowerTol = D.SteadyMotorCurrent-D.SteadyMotorCurrent*Tol;
            hold on 
            line(D.loc.Ia_x.Peak,D.loc.Ia_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                        % point for PeakCurrent
            line([D.loc.Ia_x.PeakOffset,D.loc.Ia_x.Peak],[D.loc.Ia_y.PeakOffset,D.loc.Ia_y.Peak],'Color',[0 204 0]/255);        % line for PeakCurrent
            line([0,D.loc.Ia_x.Peak],[D.loc.Ia_y.Peak,D.loc.Ia_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % horizontal line for PeakCurrent
            line([D.loc.Ia_x.Peak,D.loc.Ia_x.Peak],[0,D.loc.Ia_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % vertical line for PeakCurrent
            line([0,Time.Stop],[D.loc.Ia_y.Steady,D.loc.Ia_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');                   % horizontal line for SteadyCurrent
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(D.loc.Ia_x.PeakOffset,D.loc.Ia_y.PeakOffset,['Peak Current = ',num2str(D.PeakMotorCurrent),' ',Figures.Units.current],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);           % label for PeakCurrent
            text(D.loc.Ia_x.Settling,D.loc.Ia_y.Steady,['Steady State Current = ',num2str(D.SteadyMotorCurrent),' ',Figures.Units.current],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);    % label for horizontal line for SteadyCurrent
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(D.Time.Vector,D.Ia);
        f_figureTitle(D, Figures.Ia.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Ia.Label);
        hold off
        
        D.figureList = [D.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            D.loc.V_x.PeakOffset = D.TimeToPeak+D.TimeToPeak*.1;
            D.loc.V_x.Peak = D.TimeToPeak;
            D.loc.V_x.Settling = D.SettlingTime;
            D.loc.V_y.PeakOffset = D.PeakVelocity-D.PeakVelocity*0.03;
            D.loc.V_y.Peak = D.PeakVelocity;
            D.loc.V_y.Steady = D.SteadyVelocity;
            D.loc.V_y.Overshoot = (D.SteadyVelocity+D.PeakVelocity)/2;
            D.loc.V_y.UpperTol = D.SteadyVelocity+D.SteadyVelocity*Tol;
            D.loc.V_y.LowerTol = D.SteadyVelocity-D.SteadyVelocity*Tol;
            hold on 
            line(D.loc.V_x.Peak,D.loc.V_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                  % point for PeakVelocity
            line([D.loc.V_x.PeakOffset,D.loc.V_x.Peak],[D.loc.V_y.PeakOffset,D.loc.V_y.Peak],'Color',[0 204 0]/255);        % line for PeakVelocity
            line([0,D.loc.V_x.Peak],[D.loc.V_y.Peak,D.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');              % horizontal line for PeakVelocity
            line([D.loc.V_x.Peak,D.loc.V_x.Peak],[0,D.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');              % vertical line for PeakVelocity
            line([0,D.Time.Stop],[D.loc.V_y.Steady,D.loc.V_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');             % horizontal line for SteadyVelocity
            line([0,D.Time.Stop],[D.loc.V_y.UpperTol,D.loc.V_y.UpperTol],'Color',[204 204 0]/255,'LineStyle','--');       % horizontal line for SteadyVelocity UpperTol
            line([0,D.Time.Stop],[D.loc.V_y.LowerTol,D.loc.V_y.LowerTol],'Color',[204 204 0]/255,'LineStyle','--');       % horizontal line for SteadyVelocity LowerTol
            line([D.loc.V_x.Settling,D.loc.V_x.Settling],[0,D.Vx(max(find(abs(((D.Vx-D.SteadyVelocity)/D.SteadyVelocity))>Tol))+1)],'Color',[255 153 51]/255,'LineStyle','--');% vertical line for SettlingTime
            line([0.1,0.1],[D.loc.V_y.Steady,D.loc.V_y.Peak],'Color',[0 0 0]/255,'LineStyle',':');                    % vertical line for Overshoot
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(D.loc.V_x.PeakOffset,D.loc.V_y.PeakOffset,['Peak Velocity = ',num2str(D.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);    % label for PeakVelocity
            text(0.001,D.loc.V_y.Peak,['Peak Velocity = ',num2str(D.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                     % label for horizontal line for PeakVelocity
            text(D.loc.V_x.Peak,0,['  Time to Peak = ',num2str(D.TimeToPeak),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                        % label for vertical line for PeakVelocity
            text(0.001,D.loc.V_y.Steady,['Steady State Velocity = ',num2str(D.SteadyVelocity),' ',Figures.Units.velocity],'VerticalAlignment','top','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);            % label for horizontal line for SteadyVelocity
            text(0.001,D.loc.V_y.LowerTol-0.01,['Steady State Velocity Range: \pm ',num2str(Tol*100),' ',Figures.Units.percent],'VerticalAlignment','top','HorizontalAlignment','left','color',[204 204 0]/255,'FontSize',8);                         % label for horizontal lines for Tolerance Range
            text(D.loc.V_x.Settling,0,['  Settling Time = ',num2str(D.SettlingTime),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 153 51]/255,'FontSize',8);             % label for vertical line for SettlingTime
            text(0.05,D.loc.V_y.Overshoot,['  % Overshoot = ',num2str(D.Overshoot),' ',Figures.Units.percent],'VerticalAlignment','cap','HorizontalAlignment','left','color',[0 0 0]/255,'FontSize',8);                         % label for vertical line for Overshoot
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(D.Time.Vector,D.Vx);
        f_figureTitle(D, Figures.Vx.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label); 
        hold off
    end
end

%% Output for E1 (MATLAB)
    E1.txtOvershoot = ['1. % overshoot:',num2str(E1.Overshoot),' ',Figures.Units.percent];
    E1.txtPeakPosition = ['2. peak position (Peak Amplitude):',num2str(E1.PeakPosition),' ',Figures.Units.distance];
    E1.txtTimeToPeak = ['3. time to peak:',num2str(E1.TimeToPeak),' ',Figures.Units.time];
    E1.txtSteadyPosition = ['4. steady state position:',num2str(E1.SteadyPosition),' ',Figures.Units.distance];
    E1.txtSettlingTime = ['5. settling time:',num2str(E1.SettlingTime),' ',Figures.Units.time];
if(E1.showPerformanceInfo)
    disp(['<<<   ',E1.txtTitle,'   >>>']);
    disp('Performance values:');
    disp(E1.txtOvershoot); 
    disp(E1.txtPeakPosition); 
    disp(E1.txtTimeToPeak); 
    disp(E1.txtSteadyPosition);
    disp(E1.txtSettlingTime);
end
if(E1.plotGraphs)
    E1.figureList = [E1.figureList f_newFigure()];
    %annotations--------------------------------------------------------------------------------------------------------------------
        E1.loc.P_x.PeakOffset = E1.TimeToPeak+E1.TimeToPeak*0.3;
        E1.loc.P_x.Peak = E1.TimeToPeak;
        E1.loc.P_x.Settling = E1.SettlingTime;
        E1.loc.P_y.PeakOffset = E1.PeakPosition-E1.PeakPosition*0.03;
        E1.loc.P_y.Peak = E1.PeakPosition;
        E1.loc.P_y.Steady = E1.SteadyPosition;
        E1.loc.P_y.Overshoot = (E1.SteadyPosition+E1.PeakPosition)/2;
        E1.loc.P_y.UpperTol = E1.SteadyPosition+E1.SteadyPosition*Tol;
        E1.loc.P_y.LowerTol = E1.SteadyPosition-E1.SteadyPosition*Tol;
        hold on 
        line(E1.loc.P_x.Peak,E1.loc.P_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                    % point for PeakPosition
        line([E1.loc.P_x.PeakOffset,E1.loc.P_x.Peak],[E1.loc.P_y.PeakOffset,E1.loc.P_y.Peak],'Color',[0 204 0]/255);        % line for PeakPosition
        line([0,E1.loc.P_x.Peak],[E1.loc.P_y.Peak,E1.loc.P_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % horizontal line for PeakPosition
        line([E1.loc.P_x.Peak,E1.loc.P_x.Peak],[0,E1.loc.P_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % vertical line for PeakPosition
        line([0,E1.Time.Stop],[E1.loc.P_y.Steady,E1.loc.P_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');               % horizontal line for SteadyPosition
        line([0,E1.Time.Stop],[E1.loc.P_y.UpperTol,E1.loc.P_y.UpperTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyPosition UpperTol
        line([0,E1.Time.Stop],[E1.loc.P_y.LowerTol,E1.loc.P_y.LowerTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyPosition LowerTol
        line([E1.loc.P_x.Settling,E1.loc.P_x.Settling],[0,E1.Step.Output(max(find(abs(((E1.Step.Output-E1.SteadyPosition)/E1.SteadyPosition))>Tol))+1)],'Color',[255 153 51]/255,'LineStyle','--');% vertical line for SettlingTime
        line([0.05,0.05],[E1.loc.P_y.Steady,E1.loc.P_y.Peak],'Color',[0 0 0]/255,'LineStyle',':');                      % vertical line for Overshoot
        %------------------------------------------------------------------------------------------------------------------------------- 
        text(E1.loc.P_x.PeakOffset,E1.loc.P_y.PeakOffset,['Peak Position = ',num2str(E1.PeakPosition),' ',Figures.Units.distance],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8); % label for PeakPosition
        text(0.001,E1.loc.P_y.Peak,['Peak Position = ',num2str(E1.PeakPosition),' ',Figures.Units.distance],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                   % label for horizontal line for PeakPosition
        text(E1.loc.P_x.Peak,2,['  Time to Peak = ',num2str(E1.TimeToPeak),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                      % label for vertical line for PeakPosition
        text(0.001,E1.loc.P_y.Steady,['Steady State Position = ',num2str(E1.SteadyPosition),' ',Figures.Units.distance],'VerticalAlignment','top','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);          % label for horizontal line for SteadyPosition
        text(0.001,E1.loc.P_y.LowerTol-0.5,['Steady State Position Range: \pm ',num2str(Tol*100),' ',Figures.Units.percent],'VerticalAlignment','top','HorizontalAlignment','left','color',[204 204 0]/255,'FontSize',8);                        % label for horizontal lines for Tolerance Range
        text(E1.loc.P_x.Settling,2,['  Settling Time = ',num2str(E1.SettlingTime),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 153 51]/255,'FontSize',8);           % label for vertical line for SettlingTime
        text(0.05,E1.loc.P_y.Overshoot,['  % Overshoot = ',num2str(E1.Overshoot),' ',Figures.Units.percent],'VerticalAlignment','cap','HorizontalAlignment','left','color',[0 0 0]/255,'FontSize',8);                       % label for vertical line for Overshoot
    %-------------------------------------------------------------------------------------------------------------------------------
    plot(E1.Time.Vector,E1.Step.Output);
    f_figureTitle(E1, Figures.Xx.Title);
    xlabel(Figures.timeLabel); ylabel(Figures.Xx.Label);
    
    
    E1.figureList = [E1.figureList f_newFigure()];
    f_figureTitle(E1, Figures.rlocus.Title);
    rlocus(E1.T.Gx);
end

%% Calculate Output for E2 (SIMULINK)
if(E2.runSimulink)
%     %Supports negative value by preserving the sign.
%     E2.PeakVelocity = max(abs(E2.Vx))*sign(max(abs(E2.Vx)));
%     %Finds 'tout' of the index of the values of 'Vx' that are PeakVelocity
%     E2.TimeToPeakVelocity = E2.tout(find(E2.Vx==E2.PeakVelocity));
%     E2.SteadyVelocity = E2.Vx(numel(E2.Vx)-1);
%     %Finds 'tout' of the index of the values of 'Vx' that are within 2% (0.02) of Steady State
%     E2.SettlingTimeVelocity = E2.tout(max(find(abs(((E2.Vx-E2.SteadyVelocity)/E2.SteadyVelocity))>Tol))+1); 
%     E2.OvershootVelocity = (E2.PeakVelocity-E2.SteadyVelocity)/E2.SteadyVelocity *100;
    
    %Supports negative value by preserving the sign.
    E2.PeakPosition = max(abs(E2.Xx))*sign(max(abs(E2.Xx)));
    %Finds 'tout' of the index of the values of 'Vx' that are PeakVelocity
    E2.TimeToPeakPosition = E2.tout(find(E2.Xx==E2.PeakPosition));
    E2.SteadyPosition = E2.Xx(numel(E2.Xx)-1);
    %Finds 'tout' of the index of the values of 'Vx' that are within 2% (0.02) of Steady State
    E2.SettlingTimePosition = E2.tout(max(find(abs(((E2.Xx-E2.SteadyPosition)/E2.SteadyPosition))>Tol))+1); 
    E2.OvershootPosition = (E2.PeakPosition-E2.SteadyPosition)/E2.SteadyPosition *100;

    E2.PeakMotorCurrent = max(E2.Ia);
    E2.TimeToPeakMotorCurrent = E2.tout(find(E2.Ia==E2.PeakMotorCurrent));
    E2.SteadyMotorCurrent = E2.Ia(numel(E2.Ia)-1); % Get last value recorded
    E2.TimeToSteadyMotorCurrent = E2.tout(find(E2.Ia==E2.SteadyMotorCurrent));
    
    E2.PeakMotorTorque = max(E2.Tm);
    E2.TimeToPeakMotorTorque = E2.tout(find(E2.Tm==E2.PeakMotorTorque));
    E2.SteadyMotorTorque = E2.Tm(numel(E2.Tm)-1); % Get last value recorded
    E2.TimeToSteadyMotorTorque = E2.tout(find(E2.Tm==E2.SteadyMotorTorque));

    % Prepare Output
%     E2.txtOvershootVelocity =       ['1. % overshoot: ',num2str(E2.OvershootVelocity),' ',Figures.Units.percent];
%     E2.txtPeakVelocity =            ['2. peak Velocity (Peak Amplitude): ',num2str(E2.PeakVelocity),' ',Figures.Units.velocity];
%     E2.txtTimeToPeakVelocity =      ['3. time to peak: ',num2str(E2.TimeToPeakVelocity),' ',Figures.Units.time];
%     E2.txtSteadyVelocity =          ['4. steady state Velocity: ',num2str(E2.SteadyVelocity),' ',Figures.Units.velocity];
%     E2.txtSettlingTimeVelocity =    ['5. settling time: ',num2str(E2.SettlingTimeVelocity),' ',Figures.Units.time];
    
    E2.txtOvershootPosition =       ['1. % overshoot: ',num2str(E2.OvershootPosition),' ',Figures.Units.percent];
    E2.txtPeakPosition =            ['2. peak position (Peak Amplitude): ',num2str(E2.PeakPosition),' ',Figures.Units.distance];
    E2.txtTimeToPeakPosition =      ['3. time to peak: ',num2str(E2.TimeToPeakPosition),' ',Figures.Units.time];
    E2.txtSteadyPosition =          ['4. steady state position: ',num2str(E2.SteadyPosition),' ',Figures.Units.distance];
    E2.txtSettlingTimePosition =    ['5. settling time: ',num2str(E2.SettlingTimePosition),' ',Figures.Units.time];

    E2.txtPeakMotorCurrent =    ['6. peak motor input current: ',num2str(E2.PeakMotorCurrent),' ',Figures.Units.current];
    E2.txtSteadyMotorCurrent =  ['7. steady-state motor input current: ',num2str(E2.SteadyMotorCurrent),' ',Figures.Units.current];
    E2.txtPeakMotorTorque =     ['8. peak motor output torque: ',num2str(E2.PeakMotorTorque),' ',Figures.Units.torque];
    E2.txtSteadyMotorTorque =   ['9. steady-state motor output torque: ',num2str(E2.SteadyMotorTorque),' ',Figures.Units.torque];
end
%% Output for E2 (SIMULINK)
if(E2.runSimulink)
    if(E2.showPerformanceInfo)
        disp(['<<<   ',E2.txtTitle,'   >>>']);
%         disp('Performance values:');
%         disp(E2.txtOvershootVelocity); 
%         disp(E2.txtPeakVelocity); 
%         disp(E2.txtTimeToPeakVelocity); 
%         disp(E2.txtSteadyVelocity);
%         disp(E2.txtSettlingTimeVelocity);

        disp(E2.txtOvershootPosition); 
        disp(E2.txtPeakPosition); 
        disp(E2.txtTimeToPeakPosition); 
        disp(E2.txtSteadyPosition);
        disp(E2.txtSettlingTimePosition);

        disp(E2.txtPeakMotorCurrent);
        disp(E2.txtSteadyMotorCurrent);
        disp(E2.txtPeakMotorTorque);
        disp(E2.txtSteadyMotorTorque);
    end
    if(E2.plotGraphs)
        E2.figureList = [E2.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            E2.loc.Tm_x.PeakOffset = E2.TimeToPeakMotorTorque+E2.TimeToPeakMotorTorque*1;
            E2.loc.Tm_x.Peak = E2.TimeToPeakMotorTorque;
            E2.loc.Tm_x.Settling = E2.TimeToSteadyMotorTorque;
            E2.loc.Tm_y.PeakOffset = E2.PeakMotorTorque-E2.PeakMotorTorque*0.03;
            E2.loc.Tm_y.Peak = E2.PeakMotorTorque;
            E2.loc.Tm_y.Steady = E2.SteadyMotorTorque;
            E2.loc.Tm_y.Overshoot = (E2.SteadyMotorTorque+E2.PeakMotorTorque)/2;
            E2.loc.Tm_y.UpperTol = E2.SteadyMotorTorque+E2.SteadyMotorTorque*Tol;
            E2.loc.Tm_y.LowerTol = E2.SteadyMotorTorque-E2.SteadyMotorTorque*Tol;
            hold on 
            line(E2.loc.Tm_x.Peak,E2.loc.Tm_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                        % point for PeakTorque
            line([E2.loc.Tm_x.PeakOffset,E2.loc.Tm_x.Peak],[E2.loc.Tm_y.PeakOffset,E2.loc.Tm_y.Peak],'Color',[0 204 0]/255);        % line for PeakTorque
            line([0,E2.loc.Tm_x.Peak],[E2.loc.Tm_y.Peak,E2.loc.Tm_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % horizontal line for PeakTorque
            line([E2.loc.Tm_x.Peak,E2.loc.Tm_x.Peak],[0,E2.loc.Tm_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % vertical line for PeakTorque
            line([0,Time.Stop],[E2.loc.Tm_y.Steady,E2.loc.Tm_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');                   % horizontal line for SteadyTorque
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(E2.loc.Tm_x.PeakOffset,E2.loc.Tm_y.PeakOffset,['Peak Torque = ',num2str(E2.PeakMotorTorque),' ',Figures.Units.torque],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);           % label for PeakTorque
            text(0,E2.loc.Tm_y.Steady,['Steady State Torque = ',num2str(E2.SteadyMotorTorque),' ',Figures.Units.torque],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);    % label for horizontal line for SteadyTorque
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(E2.Time.Vector,E2.Tm);
        f_figureTitle(E2, Figures.Tm.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Tm.Label);
        hold off
        
        E2.figureList = [E2.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            E2.loc.Ia_x.PeakOffset = E2.TimeToPeakMotorCurrent+E2.TimeToPeakMotorCurrent*1;
            E2.loc.Ia_x.Peak = E2.TimeToPeakMotorCurrent;
            E2.loc.Ia_x.Settling = E2.TimeToSteadyMotorCurrent;
            E2.loc.Ia_y.PeakOffset = E2.PeakMotorCurrent-E2.PeakMotorCurrent*0.03;
            E2.loc.Ia_y.Peak = E2.PeakMotorCurrent;
            E2.loc.Ia_y.Steady = E2.SteadyMotorCurrent;
            E2.loc.Ia_y.Overshoot = (E2.SteadyMotorCurrent+E2.PeakMotorCurrent)/2;
            E2.loc.Ia_y.UpperTol = E2.SteadyMotorCurrent+E2.SteadyMotorCurrent*Tol;
            E2.loc.Ia_y.LowerTol = E2.SteadyMotorCurrent-E2.SteadyMotorCurrent*Tol;
            hold on 
            line(E2.loc.Ia_x.Peak,E2.loc.Ia_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                        % point for PeakCurrent
            line([E2.loc.Ia_x.PeakOffset,E2.loc.Ia_x.Peak],[E2.loc.Ia_y.PeakOffset,E2.loc.Ia_y.Peak],'Color',[0 204 0]/255);        % line for PeakCurrent
            line([0,E2.loc.Ia_x.Peak],[E2.loc.Ia_y.Peak,E2.loc.Ia_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % horizontal line for PeakCurrent
            line([E2.loc.Ia_x.Peak,E2.loc.Ia_x.Peak],[0,E2.loc.Ia_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');                 % vertical line for PeakCurrent
            line([0,Time.Stop],[E2.loc.Ia_y.Steady,E2.loc.Ia_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');                   % horizontal line for SteadyCurrent
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(E2.loc.Ia_x.PeakOffset,E2.loc.Ia_y.PeakOffset,['Peak Current = ',num2str(E2.PeakMotorCurrent),' ',Figures.Units.current],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);           % label for PeakCurrent
            text(0,E2.loc.Ia_y.Steady,['Steady State Current = ',num2str(E2.SteadyMotorCurrent),' ',Figures.Units.current],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);    % label for horizontal line for SteadyCurrent
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(E2.Time.Vector,E2.Ia);
        f_figureTitle(E2, Figures.Ia.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Ia.Label);
        hold off
        
%         E2.figureList = [E2.figureList f_newFigure()];
%         %annotations--------------------------------------------------------------------------------------------------------------------
%             E2.loc.V_x.PeakOffset = E2.TimeToPeakVelocity+E2.TimeToPeakVelocity*0.3;
%             E2.loc.V_x.Peak = E2.TimeToPeakVelocity;
%             E2.loc.V_x.Settling = E2.SettlingTimeVelocity;
%             E2.loc.V_y.PeakOffset = E2.PeakVelocity-E2.PeakVelocity*0.03;
%             E2.loc.V_y.Peak = E2.PeakVelocity;
%             E2.loc.V_y.Steady = E2.SteadyVelocity; 
%             E2.loc.V_y.Overshoot = (E2.SteadyVelocity+E2.PeakVelocity)/2;
%             E2.loc.V_y.UpperTol = E2.SteadyVelocity+E2.SteadyVelocity*Tol;
%             E2.loc.V_y.LowerTol = E2.SteadyVelocity-E2.SteadyVelocity*Tol;
%             hold on 
%             line(E2.loc.V_x.Peak,E2.loc.V_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                    % point for PeakVelocity
%             line([E2.loc.V_x.PeakOffset,E2.loc.V_x.Peak],[E2.loc.V_y.PeakOffset,E2.loc.V_y.Peak],'Color',[0 204 0]/255);        % line for PeakVelocity
%             line([0,E2.loc.V_x.Peak],[E2.loc.V_y.Peak,E2.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % horizontal line for PeakVelocity
%             line([E2.loc.V_x.Peak,E2.loc.V_x.Peak],[0,E2.loc.V_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % vertical line for PeakVelocity
%             line([0,E2.Time.Stop],[E2.loc.V_y.Steady,E2.loc.V_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');               % horizontal line for SteadyVelocity
%             line([0,E2.Time.Stop],[E2.loc.V_y.UpperTol,E2.loc.V_y.UpperTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyVelocity UpperTol
%             line([0,E2.Time.Stop],[E2.loc.V_y.LowerTol,E2.loc.V_y.LowerTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyVelocity LowerTol
%             line([E2.loc.V_x.Settling,E2.loc.V_x.Settling],[0,E2.Vx(max(find(abs(((E2.Vx-E2.SteadyVelocity)/E2.SteadyVelocity))>Tol))+1);],'Color',[255 153 51]/255,'LineStyle','--');% vertical line for SettlingTime
%             line([0.05,0.05],[E2.loc.V_y.Steady,E2.loc.V_y.Peak],'Color',[0 0 0]/255,'LineStyle',':');                      % vertical line for Overshoot
%             %------------------------------------------------------------------------------------------------------------------------------- 
%             text(E2.loc.V_x.PeakOffset,E2.loc.V_y.PeakOffset,['Peak Velocity = ',num2str(E2.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8); % label for PeakVelocity
%             text(0.001,E2.loc.V_y.Peak,['Peak Velocity = ',num2str(E2.PeakVelocity),' ',Figures.Units.velocity],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                   % label for horizontal line for PeakVelocity
%             text(E2.loc.V_x.Peak,100,['  Time to Peak = ',num2str(E2.TimeToPeakVelocity),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                      % label for vertical line for PeakVelocity
%             text(0.001,E2.loc.V_y.Steady,['Steady State Velocity = ',num2str(E2.SteadyVelocity),' ',Figures.Units.velocity],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);          % label for horizontal line for SteadyVelocity
%             text(0.001,E2.loc.V_y.LowerTol,['Steady State Velocity Range: \pm ',num2str(Tol*100),' ',Figures.Units.percent],'VerticalAlignment','top','HorizontalAlignment','left','color',[204 204 0]/255,'FontSize',8);                        % label for horizontal lines for Tolerance Range
%             text(E2.loc.V_x.Settling,0.1,['  Settling Time = ',num2str(E2.SettlingTimeVelocity),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 153 51]/255,'FontSize',8);           % label for vertical line for SettlingTime
%             text(0.05,E2.loc.V_y.Overshoot,['  % Overshoot = ',num2str(E2.OvershootVelocity),' ',Figures.Units.percent],'VerticalAlignment','cap','HorizontalAlignment','left','color',[0 0 0]/255,'FontSize',8);                       % label for vertical line for Overshoot
%         %-------------------------------------------------------------------------------------------------------------------------------
%         plot(E2.Time.Vector,E2.Vx);
%         f_figureTitle(E2, Figures.Vx.Title);
%         xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);
%         hold off
        
        E2.figureList = [E2.figureList f_newFigure()];
        %annotations--------------------------------------------------------------------------------------------------------------------
            E2.loc.P_x.PeakOffset = E2.TimeToPeakPosition+E2.TimeToPeakPosition*0.3;
            E2.loc.P_x.Peak = E2.TimeToPeakPosition;
            E2.loc.P_x.Settling = E2.SettlingTimePosition;
            E2.loc.P_y.PeakOffset = E2.PeakPosition-E2.PeakPosition*0.03;
            E2.loc.P_y.Peak = E2.PeakPosition;
            E2.loc.P_y.Steady = E2.SteadyPosition; 
            E2.loc.P_y.Overshoot = (E2.SteadyPosition+E2.PeakPosition)/2;
            E2.loc.P_y.UpperTol = E2.SteadyPosition+E2.SteadyPosition*Tol;
            E2.loc.P_y.LowerTol = E2.SteadyPosition-E2.SteadyPosition*Tol;
            hold on 
            line(E2.loc.P_x.Peak,E2.loc.P_y.Peak,'Color',[0 204 0]/255,'LineStyle','.','MarkerSize',15);                    % point for PeakPosition
            line([E2.loc.P_x.PeakOffset,E2.loc.P_x.Peak],[E2.loc.P_y.PeakOffset,E2.loc.P_y.Peak],'Color',[0 204 0]/255);        % line for PeakPosition
            line([0,E2.loc.P_x.Peak],[E2.loc.P_y.Peak,E2.loc.P_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % horizontal line for PeakPosition
            line([E2.loc.P_x.Peak,E2.loc.P_x.Peak],[0,E2.loc.P_y.Peak],'Color',[0 204 0]/255,'LineStyle','--');               % vertical line for PeakPosition
            line([0,E2.Time.Stop],[E2.loc.P_y.Steady,E2.loc.P_y.Steady],'Color',[255 0 0]/255,'LineStyle','--');               % horizontal line for SteadyPosition
            line([0,E2.Time.Stop],[E2.loc.P_y.UpperTol,E2.loc.P_y.UpperTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyPosition UpperTol
            line([0,E2.Time.Stop],[E2.loc.P_y.LowerTol,E2.loc.P_y.LowerTol],'Color',[204 204 0]/255,'LineStyle','--');         % horizontal line for SteadyPosition LowerTol
            line([E2.loc.P_x.Settling,E2.loc.P_x.Settling],[0,E2.Xx(max(find(abs(((E2.Xx-E2.SteadyPosition)/E2.SteadyPosition))>Tol))+1);],'Color',[255 153 51]/255,'LineStyle','--');% vertical line for SettlingTime
            line([0.05,0.05],[E2.loc.P_y.Steady,E2.loc.P_y.Peak],'Color',[0 0 0]/255,'LineStyle',':');                      % vertical line for Overshoot
            %------------------------------------------------------------------------------------------------------------------------------- 
            text(E2.loc.P_x.PeakOffset,E2.loc.P_y.PeakOffset,['Peak Position = ',num2str(E2.PeakPosition),' ',Figures.Units.distance],'VerticalAlignment','top','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8); % label for PeakPosition
            text(0.001,E2.loc.P_y.Peak,['Peak Position = ',num2str(E2.PeakPosition),' ',Figures.Units.distance],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                   % label for horizontal line for PeakPosition
            text(E2.loc.P_x.Peak,2,['  Time to Peak = ',num2str(E2.TimeToPeakPosition),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8);                      % label for vertical line for PeakPosition
            text(0.001,E2.loc.P_y.Steady,['Steady State Position = ',num2str(E2.SteadyPosition),' ',Figures.Units.distance],'VerticalAlignment','top','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8);          % label for horizontal line for SteadyPosition
            text(0.001,E2.loc.P_y.LowerTol-0.5,['Steady State Position Range: \pm ',num2str(Tol*100),' ',Figures.Units.percent],'VerticalAlignment','top','HorizontalAlignment','left','color',[204 204 0]/255,'FontSize',8);                        % label for horizontal lines for Tolerance Range
            text(E2.loc.P_x.Settling,2,['  Settling Time = ',num2str(E2.SettlingTimePosition),' ',Figures.Units.time],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 153 51]/255,'FontSize',8);           % label for vertical line for SettlingTime
            text(0.05,E2.loc.P_y.Overshoot,['  % Overshoot = ',num2str(E2.OvershootPosition),' ',Figures.Units.percent],'VerticalAlignment','cap','HorizontalAlignment','left','color',[0 0 0]/255,'FontSize',8);                       % label for vertical line for Overshoot
        %-------------------------------------------------------------------------------------------------------------------------------
        plot(E2.Time.Vector,E2.Xx);
        f_figureTitle(E2, Figures.Xx.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Xx.Label); 
        hold off
    end
end

%% Compare C3 and D1   
if( (C2.runSimulink) && (C3.plotGraphs))
    C3.Vx_Matlab = C1.Step.Output;
    C3.Vx_Simulink = C2.Vx;
    C3.MatlabMinusSimulink = C3.Vx_Matlab - C3.Vx_Simulink;
    
    C3.figureList = [C3.figureList f_newFigure()];
    plot(C3.Time.Vector,C3.Vx_Matlab,'r');
    hold on
    plot(C3.Time.Vector,C3.Vx_Simulink+0.1,'b');
    hold off
    title({'C3a. Overlay Matlab and Simulink for Open-Loop (Neglecting Gravity)',Figures.Vx.Title});
    xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);
    legend('Matlab', 'Simulink','Location','Best');
    
    C3.figureList = [C3.figureList f_newFigure()];
    plot(C3.Time.Vector,C3.MatlabMinusSimulink);
    title({'C3b. Matlab Minus Simulink for Open-Loop (Neglecting Gravity)',Figures.Vx.Title});
    xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);
end
if( (C2.runSimulink) && (D.runSimulink) && (D1.plotGraphs))
    D1.Vx_NoGravity = C2.Vx;
    D1.Vx_Gravity = D.Vx;
    D1.NoGravityMinusGravity = D1.Vx_NoGravity - D1.Vx_Gravity;
    
    D1.figureList = [D1.figureList f_newFigure()];
    plot(D1.Time.Vector,D1.Vx_Gravity,'r');
    hold on
    plot(D1.Time.Vector,D1.Vx_NoGravity,'b');
    hold off
    title({'D1a.Overlay Simulink with and without Gravity',Figures.Vx.Title});
    xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);
    legend('With Gravity', 'Without Gravity','Location','Best');
    
    D1.figureList = [D1.figureList f_newFigure()];
    plot(D1.Time.Vector,D1.NoGravityMinusGravity);
    title({'D1b. Simulink: Without-Gravity Minus With-Gravity',Figures.Vx.Title});
    xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);
end

%% Find all figures in workspace
Figures.NumberOfFigures = gcf;
% Figure Positions
    Figures.rightScreen_Small_normalized=[0.525, 0.075, 0.45, 0.825]; %[x y w h]
    Figures.leftScreen_Small_normalized=[0.025, 0.075, 0.45, 0.825]; %[x y w h]
    Figures.rightScreen_fit_normalized=[0.505, 0.05, 0.492, 0.8725]; %[x y w h]
    Figures.leftScreen_fit_normalized=[0.005, 0.05, 0.492, 0.8725]; %[x y w h]
    Figures.fullScreen_fit_normalized=[0.005, 0.05, 0.99, 0.8725]; %[x y w h]
    Figures.wholeScreen_normalized=[0,0,1,1]; %[x y w h]
for n = 1:Figures.NumberOfFigures
    %set(n,'Units','normalized','Position',Figures.wholeScreen_normalized)
    set(n,'PaperPositionMode','auto') %Required to make plots look nice
    %Figures.Plots(n).axes = findall(n,'type','axes');
    %Figures.Plots(n).title = get(get(Figures.Plots(n).axes,'Title'),'string')
    if(isempty(get(get(get(n,'currentAxes'),'Title'),'string')))
        Figures.Plots(n).title = {E1.txtTitle(), Figures.rlocus.Title}; % Root locus plot can't be found
    else
        Figures.Plots(n).title = get(get(get(n,'currentAxes'),'Title'),'string');
    end
end
%% Save Graphs
f_savePlot =@(X,h) f_savePng(h,[X.ShortTitle,'_',Figures.Plots(h).title{2}]);
if(SavePNGImages)
    %Save C1
    for i = 1:numel(C1.figureList)
        f_savePlot(C1,C1.figureList(i));
    end
    if(C2.plotGraphs)
        %Save C2
        for i = 1:numel(C2.figureList)
            f_savePlot(C2,C2.figureList(i));
        end
    end
    if(C3.plotGraphs)
        %Save C3
        for i = 1:numel(C3.figureList)
            f_savePlot(C3,C3.figureList(i));
        end
    end
    if(D.plotGraphs)
        %Save D
        for i = 1:numel(D.figureList)
            f_savePlot(D,D.figureList(i));
        end
    end
    %Save E1
    for i = 1:numel(E1.figureList)
        f_savePlot(E1,E1.figureList(i));
    end
    if(E2.plotGraphs)
        %Save E2
        for i = 1:numel(E2.figureList)
            f_savePlot(E2,E2.figureList(i));
        end
    end
end
%% Interactive User Input
commandwindow
disp('--------------------------------------------------');
disp('EE315 Linear Control Systems Project 1 Part 2')
disp('Written By : Nathan G and Patrick S')
disp('EE315_Project1_Part2 Date: 03/16/2015 ')
disp('--------------------------------------------------');
disp('---- Interactive (N to quit,# to show figure) ----');
disp('Press:"1","Enter","4","Enter","7",..."8","13","N" ');
disp('--------------------------------------------------');
for n = 1:Figures.NumberOfFigures
        fprintf('(%2d) %-30s%s\n',n,Figures.Plots(n).title{2},Figures.Plots(n).title{1}); 
end
disp('--------------------------------------------------');
User.Done = 0;
while(User.Done==0)
    User.stringInput = input('Input:','s');
    if isempty(User.stringInput) %Default
        User.stringInput = 'N';
    end
    User.InputIsNum = (all(isstrprop(User.stringInput, 'digit'))==1);
    User.InputIsChar = (all(isstrprop(User.stringInput, 'alpha'))==1);
    if(User.InputIsNum)
        figure(str2double(User.stringInput));
    elseif(User.InputIsChar)
        switch User.stringInput(1)
            case 'N'
                User.Done = 1;
        end
    end 
end

if(clean==1)
    clear f_evalAtTime H User f_evalAtTime f_figureTitle f_newFigure f_savePlot f_savePng i clean u index j ans coeff rows n delta e epsilon k Time a0 a1 a2 a3 b0 RH_Table RH_real RH_sym k_equ 
end
warning ('off','all')
