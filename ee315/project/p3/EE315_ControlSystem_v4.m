function [Output] = EE315_ControlSystem_v4(System)
% Analysis for EE315 Project at SDSU. Spring 2015.
%   Detailed explanation goes here

%{
System:
    bGravity, bStepNotRamp, bCompensator, bLimiters,bMatlabNotSimulink,
    bAnnotate
    massElevator,timeStop
    StepAmplitude, RampSlope
    saveDir, shortName, longName, choosePlots
%}
%% Function Definitions

%% Parameter Defaults
% if(~any(strcmp('bMatlabNotSimulink',fieldnames(System))))
%     System.bMatlabNotSimulink = 0;
% end
% if(~any(strcmp('bStepNotRamp',fieldnames(System))))
%     System.bStepNotRamp = 1;
% end
% if(~any(strcmp('bLimiters',fieldnames(System))))
%     System.bLimiters = 1;
% end
% if(~any(strcmp('bGravity',fieldnames(System))))
%     System.bGravity = 1;
% end
% if(~any(strcmp('bCompensator',fieldnames(System))))
%     System.bCompensator = 1;
% end
% if(~any(strcmp('bAnnotate',fieldnames(System))))
%     System.bAnnotate = 1;
% end
% if(~any(strcmp('saveDir',fieldnames(System))))
%     System.saveDir = 'EE315_ProjectOutput';
% end
% if(~any(strcmp('shortName',fieldnames(System))))
%     System.shortName = 'NoName';
% end
% if(~any(strcmp('longName',fieldnames(System))))
%     System.longName = 'Linear Controls Project';
% end
% if(~any(strcmp('RampIntial',fieldnames(System))))
%     System.RampIntial = 0;
% end

%% Globals
global Time Gc_num Gc_den K Kmd Height InputRamp_slope InputRamp_ivel ...
    SimulinkInputSelector UN Nm Ng Tg Deq J kt kb La Ra maxPosition r...
    maxVoltageOfMotorDrive maxCurrentOfMotor maxTorqueOfMotor maxMotorDriveInput ;

%% System Constants
    Out.ToleranceForTsTr =0.02;
    System.MotorNum = 12;
%% Passed in from System input
    Me = System.massElevator;
    if(System.bStepNotRamp)  
        System.In.m_slope = 0;
        System.In.b_intial =System.StepAmplitude;
    else
        System.In.m_slope = System.RampSlope;
        System.In.b_intial= System.RampInitial;
    end
    if(~System.bCompensator)
        System.Gc = tf(1,1);
    end
%% Motor Properties
if(System.MotorNum == 12)
    Jm = 11;        % Kg*m^2
    Dm = 0.117;     % N*m*s*(rad^-1)
    Ra = 54.8 *10^-3;      % Ohms
    La = 0.92 *10^-3;      % H
    UN = 470;       % Vdc
    kt = 3285/675;      % V*(rad^-1)*s      %kt=T/IN; IN = 675, T=3285
    kb = (UN-Ra*675)/(825*(2*pi/60));% N*(m^-1) %kb=(UN-Ra*IN)/(n(2pi/60)) % n=825

    if(System.bLimiters)
        maxVoltageOfMotorDrive = 1410;
        maxCurrentOfMotor = 9999999999;%4050;%4050; %Cant be limited according to Jordan
        maxTorqueOfMotor = 9999999999;
        maxMotorDriveInput = 12; 
        maxPosition = 9999999999;
        %maxPosition = 61 - 3; %Shaft Height = 61, HeightOfCabin = 3;
    else
        maxVoltageOfMotorDrive = 9999999999;
        maxCurrentOfMotor = 9999999999;
        maxTorqueOfMotor = 9999999999;
        maxMotorDriveInput = 9999999999; 
        maxPosition = 9999999999;
    end
end
%% Physical Properties
% Dr. Hietpas Variables (For all motors)
    %assignin(ws, 'var', val)
    Jg = 250;       % Kg*m^2
    r = 1;          % m
    Dg = 75.8;      % N*m*s*(rad^-1)
    De = 151.6;     % N*s*(m^-1) 
    Nm = 60;        % Gear Teeth 
    Ng = 1440;      % Gear Teeth
    Ag = 9.82;
    
    % Depend on Me
    Deq = Dm + 4*Dg*(Nm/Ng)^2+2*De*r^2*(Nm/Ng)^2;
    J = Jm + 4*Jg*(Nm/Ng)^2+Me*r^2*(Nm/Ng)^2;
%% Account for Gravity (SIMULINK ONLY)
if(System.bGravity)
    Tg= (Nm/Ng)*r*Ag*Me;
else
    Tg = 0;
end
%% Output Setup
if ~exist(System.saveDir, 'dir')
    mkdir(System.saveDir);
end
f_figureTitle = @(X,str) title({X.longName,str});
%f_figureTitle = @(X,str) title({X.longName,[X.shortName,':',str]});
f_figureTitle_sub = @(X,str) title(str);
%f_figureTitle_sub = @(X,str) title([X.shortName,':',str]);
f_figureTitle_super = @(X) suptitle(X.longName);

f_savePng = @(H,name) print(H,[System.saveDir,'\',name,'.png'],'-dpng','-r150');
f_savePlot =@(H,name) f_savePng(H,[System.shortName,'_',name]);
f_newFigure = @() figure('units','normalized','outerposition',[0.005, 0.05, 0.99, 0.8725]);
f_moveFigure = @(H) set(figure(H), 'WindowStyle', 'docked'); % Dock window in order
try
    Figures = System.Figures;
catch 
    Figures.am.showPlot = 0;
    Figures.Wg.showPlot = 0;
    Figures.Va.showPlot = 0;
    Figures.Vb.showPlot = 0;
    if(System.bMatlabNotSimulink)
        Figures.Ia.showPlot = 0;
        Figures.Tm.showPlot = 0;
        Figures.Vx.showPlot = 0;
        Figures.VaVbIaTm.showPlot = 0;
        Figures.AmWgVxXx.showPlot = 0;
        
        Figures.Xx.showPlot = 1;
        Figures.Bode.showPlot = 1;
        Figures.RLocus.showPlot = 1;
    end
    if(~System.bMatlabNotSimulink)
        Figures.Ia.showPlot = 1;
        Figures.Tm.showPlot = 1;
        Figures.Vx.showPlot = 1;
        Figures.Xx.showPlot = 1;

        Figures.Bode.showPlot = 1;
        Figures.RLocus.showPlot = 1;
        Figures.VaVbIaTm.showPlot = 1;
        Figures.AmWgVxXx.showPlot = 1;
    end
end
%% Figure Labels and Units
    Figures.timeLabel = 't [s]';
    Figures.Va.Title ='Voltage Supplied to Armature';
    Figures.Va.Label ='V_a [Volts]'; 
    Figures.Ia.Title ='Current through the armature';
    Figures.Ia.Label ='I_a [Amp]';
    Figures.Tm.Title ='Torque of the motor';
    Figures.Tm.Label ='T_m [Nm]';
    Figures.Am.Title ='Radial Accel. of the Motor';
    Figures.Am.Label ='\alpha_m [m/s^2]';
    Figures.Wg.Title ='Radial Velocity of the Gears';
    Figures.Wg.Label ='\omega_g [rad/s]';
    Figures.Vx.Title ='Velocity of Elevator';
    Figures.Vx.Label ='V_x [m/s]';
    Figures.Xx.Title ='Position of the Elevator';
    Figures.Xx.Label ='X_x [m]';
    Figures.Vb.Title ='Back EMF from PMDC motor';
    Figures.Vb.Label ='V_b [Volts]';
    Figures.RLocus.Title = 'Root Locus';
    Figures.pz_open.Title = 'Pole Zero of KGH';
    Figures.Bode.Title = 'Bode of G,KGc,KGcG';
    Figures.RLocus.Title = 'Root Locus of KGH';

    Figures.Units.distance = '[m]';
    Figures.Units.time = '[s]';
    Figures.Units.percent = '%';
    Figures.Units.velocity = '[m/s]';
    Figures.Units.torque = '[N*m]';
    Figures.Units.current = '[A]';
%%  Open Loop TF
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

    System.G = Gx.TF;
    
%% Time interval for Simulation
    Time.Start = 0;
    Time.Step = 20/1000;
    Time.Stop = 20-Time.Step;
    Time.Vector = [Time.Start : Time.Step : Time.Stop];
    Out.Time = Time;
%% CALCULATE
%% Matlab
if(System.bMatlabNotSimulink)
    System.FwdLoop = System.k * series(System.G,System.Gc);
    System.T = feedback(System.FwdLoop,tf(1,1)); 
    System.T.InputUnit = 'position';
    System.T.OutputUnit = 'position';
    
    Out.dynamicSim.Input = System.In.m_slope * Time.Vector + System.In.b_intial;
    Out.dynamicSim.Output = lsim(System.T, Out.dynamicSim.Input,Time.Vector);
    Out.Xx = Out.dynamicSim.Output;
    
    if(System.bStepNotRamp)
        Out.Step.Output = step(System.T,Time.Vector,stepDataOptions('StepAmplitude',System.StepAmplitude));
        Out.Xx = Out.Step.Output;
    end % if(System.bStepNotRamp)
end % if(System.bMatlabNotSimulink)
%% Simulink Setup
if(~System.bMatlabNotSimulink)
    %% Simulink Variables
    Simulink.ModelName='EE315_Project1_Part3_Model_v2';
    Gc_num = System.Gc.num{1}; 
    Gc_den = System.Gc.den{1}; 
    K = System.k;
    Height = System.StepAmplitude;
    InputRamp_slope = System.RampSlope;
    InputRamp_ivel = System.RampInitial;
    % Simulink Options that change by system
    if(System.bStepNotRamp)
       if(System.bCompensator)
           SimulinkInputSelector = 3;
       else
           SimulinkInputSelector = 2;
       end
    else
        if(System.bCompensator)
           SimulinkInputSelector = 4;
        end
    end
    System.SimulinkInputSelector =SimulinkInputSelector;
    %% Run Simulink
    load_system(Simulink.ModelName);
    Simulink.BlockPaths = find_system(Simulink.ModelName,'Type','Block');
    Out.SimOut = sim(Simulink.ModelName, 'ReturnWorkspaceOutputs', 'on');
    % Extract Variables from Simulink
    Out.tout = Out.SimOut.get('tout');
    Out.Va = Out.SimOut.get('Va');
    Out.Ia = Out.SimOut.get('Ia');
    Out.Tm = Out.SimOut.get('Tm');
    Out.Am = Out.SimOut.get('Am');
    Out.Wm = Out.SimOut.get('Wm');
    Out.Wg = Out.SimOut.get('Wg');
    Out.Ax = Out.SimOut.get('Ax');
    Out.Vx = Out.SimOut.get('Vx');
    Out.Xx = Out.SimOut.get('Xx');
    Out.Vb = Out.SimOut.get('Vb');
end %if(~bMatlabNotSimulink)

%% ANALYZE
if(System.bMatlabNotSimulink)
    if(System.bStepNotRamp)
        Out.tout = Time.Vector;
        Out.Step.Info = stepinfo(Out.Step.Output, Time.Vector);
        Out.Overshoot = Out.Step.Info.Overshoot;
        Out.PeakPosition = Out.Step.Info.Peak;
        Out.PeakTime = Out.Step.Info.PeakTime;
        Out.TimeToPeakPosition = Out.PeakTime;
        Out.SettlingTime = Out.Step.Info.SettlingTime; 
        Out.SettlingTimePosition = Out.SettlingTime;
        
        Out.SteadyPosition          = Out.Xx(numel(Out.Xx));
        Out.SettlingPosition        = Out.Xx(max(find(abs(((Out.Xx-Out.SteadyPosition)/Out.SteadyPosition))>Out.ToleranceForTsTr))+1);
        Out.RisePosition            = Out.Xx(min(find(abs(((Out.Xx-Out.SteadyPosition)/Out.SteadyPosition))<Out.ToleranceForTsTr))+1);
        Out.RiseTime                = Out.tout(min(find(abs(((Out.Xx-Out.SteadyPosition)/Out.SteadyPosition))<Out.ToleranceForTsTr))+1);    
    else    
    end
    %% Check how it does on specs 
    Out.Checklist.OvershootIsAbove5m = (System.DesignSpec.minOS<Out.SteadyPosition*Out.Overshoot);
    Out.Checklist.ElevatorNotExceedMaxPosition = (Out.PeakPosition< 58);
    Out.Checklist.TrNotLow = (System.DesignSpec.minTr<Out.RiseTime);
    Out.Checklist.TrNotHigh = (Out.RiseTime<System.DesignSpec.maxTr);
    Out.Checklist.TsNotTooHigh = (Out.RiseTime<System.DesignSpec.maxTs);
    Out.Checklist.SystemGood =Out.Checklist.OvershootIsAbove5m* ...
        Out.Checklist.ElevatorNotExceedMaxPosition*Out.Checklist.TrNotLow*Out.Checklist.TrNotHigh*Out.Checklist.TsNotTooHigh;
    Out.Checklist.sum =Out.Checklist.OvershootIsAbove5m+ ...
        Out.Checklist.ElevatorNotExceedMaxPosition+Out.Checklist.TrNotLow+Out.Checklist.TrNotHigh+Out.Checklist.TsNotTooHigh;
        
    %Out.Checklist.SpeedNotTooFast = max(abs(Out.Vx))<System.DesignSpec.maxSpeed;
    %Out.Checklist.AccelNotTooFast = max(abs(Out.Ax))<System.DesignSpec.maxAccel;
    %Out.Checklist.CurrentNotExceed6xNominal = (max(abs(Out.Ia)) < System.DesignSpec.maxCurrent);
%     Out.Checklist.SystemGood = Out.Checklist.OvershootIsAbove5m*Out.Checklist.ElevatorNotExceedMaxPosition*Out.Checklist.TrNotLow...
%         *Out.Checklist.TrNotHigh*Out.Checklist.TsNotTooHigh*Out.Checklist.SpeedNotTooFast*Out.Checklist.AccelNotTooFast*Out.Checklist.CurrentNotExceed6xNominal;
end %if(System.bMatlabNotSimulink)
if(~System.bMatlabNotSimulink)
    %Supports negative value by preserving the sign.
    Out.PeakPosition            = max(abs(Out.Xx))*sign(max(abs(Out.Xx)));
    %Finds 'tout' of the index of the values of 'Vx' that are PeakVelocity
    Out.PeakTime                = Out.tout(find(Out.Xx==Out.PeakPosition,1));
    Out.SteadyPosition          = Out.Xx(numel(Out.Xx));
    Out.SettlingPosition        = Out.Xx(max(find(abs(((Out.Xx-Out.SteadyPosition)/Out.SteadyPosition))>Out.ToleranceForTsTr))+1);
    %Finds 'tout' of the index of the values of 'Vx' that are within 2% (0.02) of Steady State
    Out.SettlingTime            = Out.tout(max(find(abs(((Out.Xx-Out.SteadyPosition)/Out.SteadyPosition))>Out.ToleranceForTsTr))+1);
    
    Out.RisePosition            = Out.Xx(min(find(abs(((Out.Xx-Out.SteadyPosition)/Out.SteadyPosition))<Out.ToleranceForTsTr))+1);
    Out.RiseTime                = Out.tout(min(find(abs(((Out.Xx-Out.SteadyPosition)/Out.SteadyPosition))<Out.ToleranceForTsTr))+1);
       
    Out.OvershootPosition       = (Out.PeakPosition-Out.SteadyPosition)/Out.SteadyPosition *100;
    Out.Overshoot               = Out.OvershootPosition;
    Out.PeakMotorCurrent        = max(Out.Ia);
    Out.PeakTimeMotorCurrent    = Out.tout(find(Out.Ia==Out.PeakMotorCurrent,1));
    Out.SteadyMotorCurrent      = Out.Ia(numel(Out.Ia)-1); % Get last value recorded
    Out.SteadyTimeMotorCurrent  = Out.tout(find(Out.Ia==Out.SteadyMotorCurrent,1));
    
    Out.PeakMotorTorque         = max(Out.Tm);
    Out.PeakTimeMotorTorque     = Out.tout(find(Out.Tm==Out.PeakMotorTorque,1));
    Out.SteadyMotorTorque       = Out.Tm(numel(Out.Tm)-1); % Get last value recorded
    Out.SteadyTimeMotorTorque   = Out.tout(find(Out.Tm==Out.SteadyMotorTorque,1));
    
    Out.txtOvershootPosition        = ['1. % overshoot: ',num2str(Out.OvershootPosition),               ' ',Figures.Units.percent];
    Out.txtPeakPosition             = ['2. peak position (Peak Amplitude): ',num2str(Out.PeakPosition), ' ',Figures.Units.distance];
    Out.txtPeakTimePosition         = ['3. time to peak: ',num2str(Out.PeakTime),                       ' ',Figures.Units.time];
    Out.txtSteadyPosition           = ['4. steady state position: ',num2str(Out.SteadyPosition),        ' ',Figures.Units.distance];
    Out.txtSettlingTimePosition     = ['5. settling time: ',num2str(Out.SettlingTime),                  ' ',Figures.Units.time];

    Out.txtPeakMotorCurrent     = ['6. peak motor input current: ',num2str(Out.PeakMotorCurrent),           ' ',Figures.Units.current];
    Out.txtSteadyMotorCurrent   = ['7. steady-state motor input current: ',num2str(Out.SteadyMotorCurrent), ' ',Figures.Units.current];
    Out.txtPeakMotorTorque      = ['8. peak motor output torque: ',num2str(Out.PeakMotorTorque),            ' ',Figures.Units.torque];
    Out.txtSteadyMotorTorque    = ['9. steady-state motor output torque: ',num2str(Out.SteadyMotorTorque),  ' ',Figures.Units.torque];
    
    %% Check how it does on specs 
    Out.Checklist.OvershootIsAbove5m = (System.DesignSpec.minOS<Out.Overshoot);
    Out.Checklist.ElevatorNotExceedMaxPosition = (Out.PeakPosition< 58);
    Out.Checklist.TrNotLow = (System.DesignSpec.minTr<Out.RiseTime);
    Out.Checklist.TrNotHigh = (Out.RiseTime<System.DesignSpec.maxTr);
    Out.Checklist.TsNotTooHigh = (Out.RiseTime<System.DesignSpec.maxTs);
    Out.Checklist.SpeedNotTooFast = max(abs(Out.Vx))<System.DesignSpec.maxSpeed;
    Out.Checklist.AccelNotTooFast = max(abs(Out.Ax))<System.DesignSpec.maxAccel;
    Out.Checklist.CurrentNotExceed6xNominal = (max(abs(Out.Ia)) < System.DesignSpec.maxCurrent);
    
    Out.Checklist.SystemGood = Out.Checklist.OvershootIsAbove5m*Out.Checklist.ElevatorNotExceedMaxPosition*Out.Checklist.TrNotLow...
        *Out.Checklist.TrNotHigh*Out.Checklist.TsNotTooHigh*Out.Checklist.SpeedNotTooFast*Out.Checklist.AccelNotTooFast*Out.Checklist.CurrentNotExceed6xNominal;
    Out.Checklist.sum = Out.Checklist.OvershootIsAbove5m+Out.Checklist.ElevatorNotExceedMaxPosition+Out.Checklist.TrNotLow...
        +Out.Checklist.TrNotHigh+Out.Checklist.TsNotTooHigh+Out.Checklist.SpeedNotTooFast+Out.Checklist.AccelNotTooFast+Out.Checklist.CurrentNotExceed6xNominal;
end %if(~System.bMatlabNotSimulink)
%% OUTPUT Plots
Figures.figureList = [];
%% Plot Tm
if(Figures.Tm.showPlot)
    Figures.figureList = [Figures.figureList f_newFigure()];
    if(System.bAnnotate==1)
        EE315_ControlSystem_Annotation_v2(Out,Time,2);
    end
    plot(Time.Vector,Out.Tm);
    f_figureTitle(System, Figures.Tm.Title);
    xlabel(Figures.timeLabel); ylabel(Figures.Tm.Label);
    hold off
    
    f_savePlot(gcf,Figures.Tm.Title);
    f_moveFigure(gcf);
end % if(Figures.Tm.showPlot)

%% Plot Ia
if(Figures.Ia.showPlot)
    Figures.figureList = [Figures.figureList f_newFigure()];
    if(System.bAnnotate==1)
        EE315_ControlSystem_Annotation_v2(Out,Time,3);
    end
    plot(Time.Vector,Out.Ia);
    f_figureTitle(System, Figures.Ia.Title);
    xlabel(Figures.timeLabel); ylabel(Figures.Ia.Label);
    hold off
    
    f_savePlot(gcf,Figures.Tm.Title);
    f_moveFigure(gcf);
end % if(Figures.Ia.showPlot)

%% Plot Xx
if(Figures.Xx.showPlot)
    Figures.figureList = [Figures.figureList f_newFigure()];
    if(System.bAnnotate==1)
        EE315_ControlSystem_Annotation_v2(Out,Time,1);
    end
    plot(Time.Vector,Out.Xx);
    f_figureTitle(System, Figures.Xx.Title);
    xlabel(Figures.timeLabel); ylabel(Figures.Xx.Label); 
    hold off
    
    f_savePlot(gcf,Figures.Xx.Title);
    f_moveFigure(gcf);
end % if(Figures.Xx.showPlot)

%% Plot Bode
if(Figures.Bode.showPlot)
    Figures.figureList = [Figures.figureList f_newFigure()];
    
    f_figureTitle(System, Figures.Bode.Title);
    %semilogx([10 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000  5000 6000 7000 8000 9000 10000 ],[2 2.3 13.3 26.3 -6.46 -11.7 -15.6 -18.8 -20.9 -23.4 -25.2 -37.9 -44.8 -50 -53.9 -57 -59.6 -62 -64 -65 ]); hold on;
    bode(System.G, System.k*System.Gc, System.k*System.Gc*System.G);
    
    legend('G','kGc','kGcG');
    
    f_savePlot(gcf,Figures.Bode.Title);
    f_moveFigure(gcf);
end % if(Figures.Bode.showPlot)

%% Plot RLocus
if(Figures.RLocus.showPlot)
    Figures.figureList = [Figures.figureList f_newFigure()];
    
    rlocus(System.k*System.Gc*System.G);
    f_figureTitle(System, Figures.RLocus.Title);
    
    f_savePlot(gcf,Figures.RLocus.Title);
    f_moveFigure(gcf);
end % if(Figures.RLocus.showPlot)

%% Plot VaVbIaTm
if(Figures.VaVbIaTm.showPlot)
    Figures.figureList = [Figures.figureList f_newFigure()];

    subplot(2,2,1);
        plot(Time.Vector,Out.Va);
        f_figureTitle_sub(System, Figures.Va.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Va.Label);
    subplot(2,2,2);
        plot(Time.Vector,Out.Vb);
        f_figureTitle_sub(System, Figures.Vb.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Vb.Label);
    subplot(2,2,3);
        plot(Time.Vector,Out.Ia);
        f_figureTitle_sub(System, Figures.Ia.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Ia.Label);
    subplot(2,2,4);
        plot(Time.Vector,Out.Tm);
        f_figureTitle_sub(System, Figures.Tm.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Tm.Label);
    f_figureTitle_super(System);
    
    f_savePlot(gcf,'MultiPlot Va+Vb+Ia+Xx');
    f_moveFigure(gcf);
end % VaVbIaTm
%% Plot AmWgVxXx
if(Figures.AmWgVxXx.showPlot)
    Figures.figureList = [Figures.figureList f_newFigure()];
    subplot(2,2,1);
        plot(Time.Vector,Out.Am);
        f_figureTitle_sub(System, Figures.Am.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Am.Label);
    subplot(2,2,2);
        plot(Time.Vector,Out.Wg);
        f_figureTitle_sub(System, Figures.Wg.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Wg.Label);
    subplot(2,2,3);
        plot(Time.Vector,Out.Vx);
        f_figureTitle_sub(System, Figures.Vx.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Vx.Label);
    subplot(2,2,4);
        plot(Time.Vector,Out.Xx);
        f_figureTitle_sub(System, Figures.Xx.Title);
        xlabel(Figures.timeLabel); ylabel(Figures.Xx.Label);
    f_figureTitle_super(System);
        
    f_savePlot(gcf,'MultiPlot Am+Wg+Vx+Xx');
    f_moveFigure(gcf);
end % AmWgVxXx


%% Prepare for Return
if(System.bMatlabNotSimulink)
    Out.txtAnalysisSummary = sprintf('%d/4 specs, OS=%g %%, Tr=%g, Ts=%g',Out.Checklist.sum,Out.Overshoot, Out.RiseTime, Out.SettlingTime);
%     Out.Checklist.sum =Out.Checklist.OvershootIsAbove5m+ ...
%         Out.Checklist.ElevatorNotExceedMaxPosition+Out.Checklist.TrNotLow+Out.Checklist.TrNotHigh+Out.Checklist.TsNotTooHigh;       
else
    Out.txtAnalysisSummary = sprintf('%d/8 specs, OS=%g %%, Tr=%g, Ts=%g, Max{Xx=%g,Vx=%g,Ax=%g,Ia=%g}',Out.Checklist.sum,Out.Overshoot, Out.RiseTime, Out.SettlingTime,max(Out.Xx),max(abs(Out.Vx)),max(abs(Out.Ax)),max(abs(Out.Ia)) );
%     Out.Checklist.sum = Out.Checklist.OvershootIsAbove5m+Out.Checklist.ElevatorNotExceedMaxPosition+Out.Checklist.TrNotLow...
%         +Out.Checklist.TrNotHigh+Out.Checklist.TsNotTooHigh+Out.Checklist.SpeedNotTooFast+Out.Checklist.AccelNotTooFast+Out.Checklist.CurrentNotExceed6xNominal;
end
Out.Report = sprintf('%s:\t%s\n\t%s\n',System.shortName, System.longName,Out.txtAnalysisSummary);
disp(Out.Report);
System.FigureList = Figures.figureList;
Output.Out = Out;
Output.System = System;
end %End function


