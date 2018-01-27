clear all;
clc;
close all;

global Time Gc_num Gc_den K Kmd Height InputRamp_slope InputRamp_ivel ...
    SimulinkInputSelector UN Nm Ng Tg Deq J kt kb La Ra maxPosition r...
    maxVoltageOfMotorDrive maxCurrentOfMotor maxTorqueOfMotor maxMotorDriveInput ;

%% Constants
Mpeople = 62 * 10;
Me_empty = 7500;
Me_half = (Me_empty+Mpeople)/2;
Me_full = Me_empty+Mpeople;

SimulationLevel = 200; % The higher the number the more simulations are run Max is 100,000
iSave = 1;

%% Functions
f_Gc = @(Pc, a) tf([1/(a*Pc) 1], [1/Pc 1]);
f_nameLagGc = @(S) sprintf('k=%.2f;Pc=%.3f;a=%.2f',S.k,S.Pc,S.a);
f_Simulate = @(S) EE315_ControlSystem_v4(S);
f_Model_plotXx =@(Model) @(M) structfun(@(x) plot(x.Out.tout, x.Out.Xx),M,'UniformOutput',false); 
f_Model_getlongName = @(x) structfun(@(x) x.System.longName,'UniformOutput',false);
f_Model_getOut = @(M) structfun(@(x) x.Out,M,'UniformOutput',false);
f_Model_getOvershoot = @(M) structfun(@(x) x.Out.Overshoot,M,'UniformOutput',true);
f_Model_getChecklist = @(M) structfun(@(x) struct2cell(x.Out.Checklist),M,'UniformOutput',false);

Analyze.Auto = 1;
Analyze.f_getAnalysisSummary =@(M) M.LagGc_full_sim.Out.txtAnalysisSummary ;
Analyze.f_getOvershoot=@(M) M.LagGc_full_sim.Out.Overshoot;
Analyze.f_getTr =@(M) M.LagGc_full_sim.Out.RiseTime;
Analyze.f_getMaxVx =@(M) max(M.LagGc_full_sim.Out.Vx);
Analyze.HalfSim.f_getAnalysisSummary =@(M) M.LagGc_half_sim.Out.txtAnalysisSummary ;
Analyze.HalfSim.f_getOvershoot=@(M) M.LagGc_half_sim.Out.Overshoot;
Analyze.HalfSim.f_getTr =@(M) M.LagGc_half_sim.Out.RiseTime;
Analyze.HalfSim.f_getMaxVx =@(M) max(M.LagGc_half_sim.Out.Vx);
Analyze.FullSim.f_getAnalysisSummary =@(M) M.LagGc_full_sim.Out.txtAnalysisSummary ;
Analyze.FullSim.f_getOvershoot=@(M) M.LagGc_full_sim.Out.Overshoot;
Analyze.FullSim.f_getTr =@(M) M.LagGc_full_sim.Out.RiseTime;
Analyze.FullSim.f_getMaxVx =@(M) max(M.LagGc_full_sim.Out.Vx);

%% Design Specs
    System.DesignSpec.OS = 12;
    System.DesignSpec.minOS = 5/50*100;
    System.DesignSpec.minTr = 2.597;
    System.DesignSpec.maxTr = 4.9;
    System.DesignSpec.maxTs = 20;
    System.DesignSpec.maxSpeed = 20;
    System.DesignSpec.maxAccel = 7.7;
    System.DesignSpec.maxCurrent = 4050;

%% Variables able to be used
    Gains.k1 = 3.5947e+03;
    Gains.k2 = 375;
    Gains.k3 = 250;
    
    StartingLagDesign.Pc = 0.05;
    StartingLagDesign.a = 4.6;
    StartingLagDesign.k = 225;
    
    HandLagDesign.k = 1206;
    HandLagDesign.Pc = 0.0072;
    HandLagDesign.a = 1.4;
    
    JordansLagDesign.k = 209.2;
    JordansLagDesign.Pc = 1.5;
    JordansLagDesign.a = 50;
%% Default System System
System.saveDir = 'EE315_Part3Results'; 
System.bAnnotate = 1;

%These are the defaults. Order Matter down below!
System.longName = 'Default System:Step(50), Matlab, Me_empty';
System.shortName = 'Default';
System.bStepNotRamp = 1;
    System.StepAmplitude = 50; 
    System.RampSlope = 50/4; 
    System.RampInitial =0;
System.bCompensator = 0;
System.bMatlabNotSimulink = 0;
    System.bGravity = 0;
    System.bLimiters = 0;
System.massElevator = Me_empty;
    System.k = Gains.k1;
%Choose Figures
    Figures.Ia.showPlot = 0;
    Figures.Tm.showPlot = 0;
    Figures.Vx.showPlot = 0;
    Figures.Xx.showPlot = 0;
    Figures.Bode.showPlot = 0;
    Figures.RLocus.showPlot = 0;
    Figures.VaVbIaTm.showPlot = 0;
    Figures.AmWgVxXx.showPlot = 0;
    System.Figures = Figures; %Comment out line to do default
DefaultResult = f_Simulate(System);
    System = rmfield(System,'Figures');

if(400000<=SimulationLevel)
    %% A1. No Gc, Me_empty, No Gravity, StepInput [NoGcNoGrav_empty](400000)
    E1 = System; % Load Defaults
    E1.longName = 'E1.Closed-Loop Step Response neglecting gravity (MATLAB)'; E1.shortName = 'E1';
        E1.bStepNotRamp = 1;     E1.bCompensator = 0; 
        E1.bMatlabNotSimulink =1;   E1.bGravity = 0;  E1.bLimiters = 0;
        E1.massElevator = Me_empty;
        E1.k = Gains.k1;
    % NoGcNoGrav-empty
        inMatlab = E1; 
        inMatlab.k = Gains.k1;
        inMatlab.shortName = 'NoGcNoGrav-empty-Mat';
        inMatlab.longName = 'Closed-Loop Step:No Gc, No Gravity, No Load (MATLAB)';
        inSimulink = E1; 
        inSimulink.k = Gains.k1;
        inSimulink.shortName = 'NoGcNoGrav-empty-Sim';
        inSimulink.longName = 'Closed-Loop Step:No Gc, No Gravity, No Load (SIMULINK)';
        inSimulink.bMatlabNotSimulink = 0;
        if(400000<SimulationLevel)
            E1 = f_Simulate(E1);
        end
        if(400000 < SimulationLevel)
            NoGcNoGrav_empty.inMatlab = f_Simulate(inMatlab);
            NoGcNoGrav_empty.inSimulink = f_Simulate(inSimulink);
        end
    %% A2. Add People Weight. No Gc, Me_full, No Gravity, StepInput [NoGcNoGrav_full](400000)
        inMatlab.k    = Gains.k1;
            inMatlab.massElevator = Me_full;
            inMatlab.shortName = 'NoGcNoGrav-full-Mat';
            inMatlab.longName = 'Closed-Loop Step:No Gc, No Gravity, Full Load (MATLAB)';
            if(400000 < SimulationLevel)
                NoGcNoGrav_full.inMatlab_k1 = f_Simulate(inMatlab);
            end
        inSimulink.massElevator = Me_full;
            inSimulink.shortName = 'NoGcNoGrav-full-Sim';
            inSimulink.longName = 'Closed-Loop Step:No Gc, No Gravity, Full Load (SIMULINK)';  
            if(400000 < SimulationLevel)
                NoGcNoGrav_full.inSimulink = f_Simulate(inSimulink);  
            end

    %% A2a. Show Effect of Limiters and Gravity. No Gc, Me_full, No Gravity, StepInput     [NoGcNoGrav_full](3000)
        inSimulink.bGravity =0; inSimulink.bLimiters =1;
        inSimulink.k = Gains.k1; 
            inSimulink.shortName = 'NoGcNoGrav-full-Sim';
            inSimulink.longName = 'Closed-Loop Step:No Gc, No Gravity, Full Load +Limiter (SIMULINK)';
            if(3000 < SimulationLevel)
                NoGcNoGrav_full.inSimulink_Lim = f_Simulate(inSimulink);
            end

        inSimulink.bGravity =1; inSimulink.bLimiters =1;
            inSimulink.shortName = 'NoGc-full-Sim+Grav';
            inSimulink.longName = 'Closed-Loop Step:No Gc, Full Load +Limiter+Gravity (SIMULINK)';
            if(3000 < SimulationLevel)
                NoGcNoGrav_full.inSimulink_k1_Lim_Grav = f_Simulate(inSimulink);
            end
        inSimulink.k = Gains.k2; 
            inSimulink.shortName = 'NoGc-full-Sim+Grav';
            if(3000 < SimulationLevel)
                NoGcNoGrav_full.inSimulink_k2_Lim_Grav = f_Simulate(inSimulink);
            end

        inSimulink.bGravity =1; inSimulink.bLimiters =0;
            inSimulink.shortName = 'NoGc-full-Sim+Grav';
            inSimulink.longName = 'Closed-Loop Step:No Gc, Full Load +Gravity (SIMULINK)';
            if(3000 < SimulationLevel)
                NoGcNoGrav_full.inSimulink_Grav = f_Simulate(inSimulink);
            end

    %% C1. Design Compensator for Matlab. GcLag, Me_full, No Gravity, StepInput [MatlabGc.LagDesigns](550)
    %     inMatlab.bCompensator = 1;
    %     inMatlab.bGravity = 0;
    %     
    %     for iDesign = 1:length(LagDesigns)
    %         MatlabGc.LagDesigns(iDesign).Pc = MatlabGc.LagDesigns(2).Pc;
    %         MatlabGc.LagDesigns(iDesign).a = MatlabGc.LagDesigns(2).a;
    %     end
    %         MatlabGc.LagDesigns(iDesign).k = Gain(iDesign)
    %     
    %     for iDesign = 1:length(LagDesigns)
    %         inMatlab.shortName = ['LagGcNoGrav',num2str(iDesign),'-full-Mat'];
    %         inMatlab.LagDesign = MatlabGc.LagDesigns(iDesign);
    %         inMatlab.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(inMatlab.LagDesign),'] No Gravity,Full Load (MATLAB)'];
    %         inMatlab.k = MatlabGc.LagDesigns(iDesign).k;
    %         inMatlab.Gc = f_Gc(inMatlab.LagDesign.Pc,inMatlab.LagDesign.a);
    %         if(550 < SimulationLevel)
    %             Models.Lag.inMatlab(iDesign) = f_Simulate(inMatlab);
    %         end
    %     end   
    %% C1a. Show Effect of Limiters and Gravity. GcLag, Me_full, StepInput [MatlabGc.LagDesigns.inSimulink](800)
        %Redo all LagDesigns from above
    %     inSimulink.bCompensator = 1;
    %     inSimulink.bLimiter = 1; 
    %     inSimulink.bGravity = 1; 
    %     for iDesign = 1:6
    %         inSimulink.shortName = ['LagGc',num2str(iDesign),'-full-Sim'];
    %         inSimulink.LagDesign = MatlabGc.LagDesigns(iDesign);
    %         inSimulink.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(inSimulink.LagDesign),'] +Limiters+Gravity,Full Load (SIMULINK)'];
    %         inSimulink.k = inSimulink.LagDesign.k;
    %         inSimulink.Gc = f_Gc(inSimulink.LagDesign.Pc,inSimulink.LagDesign.a);
    %         if(1000 < SimulationLevel)
    %             MatlabGc.inSimulink(iDesign) = f_Simulate(inSimulink);
    %         end
    %     end
    %% D1. Design Compensator for Simulink. Me_full, With Gravity, StepInput
    %     inSimulink.bCompensator = 1;
    %     inSimulink.bGravity = 1;
    %     inSimulink.bLimiter = 1; 
    % 
    %     SimulinkGc.LagDesigns(1).k = 1206;
    %     SimulinkGc.LagDesigns(1).Pc = 0.0072;
    %     SimulinkGc.LagDesigns(1).a = 1.4;
    %     
    %     SimulinkGc.LagDesigns(2).k = 225;
    %     SimulinkGc.LagDesigns(2).Pc = 0.05;
    %     SimulinkGc.LagDesigns(2).a = 4.6;
    %     
    %     for iDesign = 3:6
    %         SimulinkGc.LagDesigns(iDesign).Pc = SimulinkGc.LagDesigns(2).Pc;
    %         SimulinkGc.LagDesigns(iDesign).a = SimulinkGc.LagDesigns(2).a;
    %     end
    %     SimulinkGc.LagDesigns(3).k = 300;
    %     SimulinkGc.LagDesigns(4).k = 325;
    %     SimulinkGc.LagDesigns(5).k = 350;
    %     SimulinkGc.LagDesigns(6).k = 375;
    % 
    %         for iDesign = 1:6
    %             inSimulink.shortName = ['LagGc',num2str(iDesign),'-full-Sim'];
    %             inSimulink.LagDesign = SimulinkGc.LagDesigns(iDesign);
    %             inSimulink.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(inSimulink.LagDesign),'] Full Load(SIMULINK)'];
    %             inSimulink.k = SimulinkGc.LagDesigns(iDesign).k;
    %             inSimulink.Gc = f_Gc(inSimulink.LagDesign.Pc,inSimulink.LagDesign.a);
    %             if(1000 < SimulationLevel)
    %             SimulinkGc.inSimulink(iDesign) = f_Simulate(inSimulink);
    %             elseif(10 < SimulationLevel)
    %                 if(iDesign ==2)
    %                    SimulinkGc.inSimulink(iDesign) = f_Simulate(inSimulink); 
    %                 end
    %             end
    %         end 
end        
%% F. Best Systems at Empty,Half, and Full Weight. Compare Ideal vs Non Ideal. (200)
if(SimulationLevel==200)
    System.bStepNotRamp = 1;
        System.StepAmplitude = 50; 
        System.RampInitial =0;
    System.bCompensator = 1;
    System.bMatlabNotSimulink = 1;
        System.bGravity = 0;
        System.bLimiters = 0;
    inMatlab = System;

    System.bMatlabNotSimulink = 0;
        System.bGravity = 1;
        System.bLimiters = 1;
        inSimulink =System;

        FinalLagDesign(1).k = 190; FinalLagDesign(1).Pc = 0.05; FinalLagDesign(1).a = 3.894; 
        FinalLagDesign(2).k = 200; FinalLagDesign(2).Pc = 0.05; FinalLagDesign(2).a = 3.6; 
        FinalLagDesign(3).k = 250; FinalLagDesign(3).Pc = 0.05; FinalLagDesign(3).a = 3.6; 
        FinalLagDesign(4).k = 300; FinalLagDesign(4).Pc = 0.05; FinalLagDesign(4).a = 3.6;
        FinalLagDesign(5).k = 350; FinalLagDesign(5).Pc = 0.05; FinalLagDesign(5).a = 3.6;
        FinalLagDesign(6).k = 400; FinalLagDesign(6).Pc = 0.05; FinalLagDesign(6).a = 3.7;
        FinalLagDesign(7).k = 1000; FinalLagDesign(7).Pc = 0.05; FinalLagDesign(7).a = 3.8;
        FinalLagDesign(8).k = 2000; FinalLagDesign(8).Pc = 0.05; FinalLagDesign(8).a = 3.9;
% FinalLagDesign(1).k = 225; FinalLagDesign(1).Pc =0.05;FinalLagDesign(1).a = 3.602; Meet Spec in FullSimulink with Current Limiter
% FinalLagDesign(1).k = 200; FinalLagDesign(1).Pc = 0.05; FinalLagDesign(1).a = 3.65; Meet Spec in HalfSimulink with Current Limiter
% FinalLagDesign(1).k = 190; FinalLagDesign(1).Pc = 0.05; FinalLagDesign(1).a = 3.88; Meet Spec in HalfSimulink without Current Limiter
% FinalLagDesign(1).k = 190; FinalLagDesign(1).Pc = 0.05; FinalLagDesign(1).a = 3.894;
    NumOfModels = 1;
    Systems=repmat(System,1,NumOfModels);
    %Models=repmat(DefaultResult,1,NumOfModels);
    for i = 1:NumOfModels  
        inMatlab.saveDir = ['EE315_Part3Results',num2str(i)];
        inSimulink.saveDir = ['EE315_Part3Results',num2str(i)];
        %FinalLagDesign(i).a = 3.6;
        %FinalLagDesign(i).Pc = 0.05;
        %FinalLagDesign(i).k = 225;
        if((i~=1) && Analyze.Auto)
            FinalLagDesign(i) = FinalLagDesign(i-1);
            %Design for 12 % OS
            if(Analyze.f_getOvershoot(Models(i-1))<System.DesignSpec.OS)
                %if Overshoot was low. Increase a
                FinalLagDesign(i).a = FinalLagDesign(i-1).a+0.001;
            else
                FinalLagDesign(i).a = FinalLagDesign(i-1).a-0.001;
            end
            
            if(Analyze.f_getTr(Models(i-1)) < System.DesignSpec.minTr)% if Tr was too high. Increase K. Otherwise decrease
                FinalLagDesign(i).k = FinalLagDesign(i-1).k+10; %Increase K if Tr too high
            elseif(System.DesignSpec.maxTr<Analyze.f_getTr(Models(i-1)))
                FinalLagDesign(i).k = FinalLagDesign(i-1).k-10; %Decrease K if Tr too low
            elseif(System.DesignSpec.maxSpeed<Analyze.f_getMaxVx(Models(i-1)))
                FinalLagDesign(i).k = FinalLagDesign(i-1).k-10; %Decrease K if speed too high
            end
            
        end
%         FSim.f_getAnalysisSummary =@(M) M.LagGc_full_sim.Out.txtAnalysisSummary ;
%         FSim.f_getOvershoot=@(M) M.LagGc_full_sim.Out.Overshoot;
%         FSim.f_getTr =@(M) M.LagGc_full_sim.Out.RiseTime;
% 
%         %% Design Specs
%             System.DesignSpec.OS = 12;
%             System.DesignSpec.minOS = 5/50*100;
%             System.DesignSpec.maxTr = 2.597;
%             System.DesignSpec.minTr = 4.9;
%             System.DesignSpec.maxTs = 20;
%             System.DesignSpec.maxSpeed = 20;
%             System.DesignSpec.maxAccel = 7.7;
%             System.DesignSpec.maxCurrent = 4050;
        
        inMatlab.LagDesign = FinalLagDesign(i);
        inMatlab.k = inMatlab.LagDesign.k;
        inMatlab.Gc = f_Gc(inMatlab.LagDesign.Pc,inMatlab.LagDesign.a);

        inSimulink.LagDesign = FinalLagDesign(i);
        inSimulink.k = inSimulink.LagDesign.k;
        inSimulink.Gc = f_Gc(inSimulink.LagDesign.Pc,inSimulink.LagDesign.a);

        Systems(i).LagGc_empty_mat = inMatlab;
            Systems(i).LagGc_empty_mat.shortName = ['LagGc-empty-Mat'];
            Systems(i).LagGc_empty_mat.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(Systems(i).LagGc_empty_mat.LagDesign),'] Empty Load(MATLAB)'];
            Systems(i).LagGc_empty_mat.massElevator = Me_empty;
        Systems(i).LagGc_half_mat = inMatlab;
            Systems(i).LagGc_half_mat.shortName = ['LagGc-half-Mat'];
            Systems(i).LagGc_half_mat.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(Systems(i).LagGc_half_mat.LagDesign),'] Half Load(MATLAB)'];
            Systems(i).LagGc_half_mat.massElevator = Me_half;
        Systems(i).LagGc_full_mat = inMatlab;
            Systems(i).LagGc_full_mat.shortName = ['LagGc-full-Mat'];
            Systems(i).LagGc_full_mat.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(Systems(i).LagGc_full_mat.LagDesign),'] Full Load(MATLAB)'];
            Systems(i).LagGc_full_mat.massElevator = Me_full;
        Systems(i).NoGc_half_mat = inMatlab;
            Systems(i).NoGc_half_mat.shortName = ['NoGc-half-Mat'];
            Systems(i).NoGc_half_mat.longName = ['Closed-Loop Step:NoGc Half Load(MATLAB)'];
            Systems(i).NoGc_half_mat.massElevator = Me_half;
            Systems(i).NoGc_half_mat.bCompensator = 0;
        Systems(i).NoGc_full_mat = inMatlab;
            Systems(i).NoGc_full_mat.shortName = ['NoGc-full-Mat'];
            Systems(i).NoGc_full_mat.longName = ['Closed-Loop Step:NoGc Full Load(MATLAB)'];
            Systems(i).NoGc_full_mat.massElevator = Me_full;
            Systems(i).NoGc_full_mat.bCompensator = 0;


        Systems(i).LagGc_empty_sim = inSimulink;
            Systems(i).LagGc_empty_sim.shortName = ['LagGc-empty-Sim'];
            Systems(i).LagGc_empty_sim.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(Systems(i).LagGc_empty_sim.LagDesign),'] Empty Load(SIMULINK)'];
            Systems(i).LagGc_empty_sim.massElevator = Me_empty;
        Systems(i).LagGc_half_sim = inSimulink;
            Systems(i).LagGc_half_sim.shortName = ['LagGc-half-Sim'];
            Systems(i).LagGc_half_sim.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(Systems(i).LagGc_half_sim.LagDesign),'] Half Load(SIMULINK)'];
            Systems(i).LagGc_half_sim.massElevator = Me_half;
        Systems(i).LagGc_full_sim = inSimulink;
            Systems(i).LagGc_full_sim.shortName = ['LagGc-full-Sim'];
            Systems(i).LagGc_full_sim.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(Systems(i).LagGc_full_sim.LagDesign),'] Full Load(SIMULINK)'];
            Systems(i).LagGc_full_sim.massElevator = Me_full;
        Systems(i).NoGc_half_sim = inSimulink;
            Systems(i).NoGc_half_sim.shortName = ['NoGc-half-Sim'];
            Systems(i).NoGc_half_sim.longName = ['Closed-Loop Step:No Gc Half Load(SIMULINK)'];
            Systems(i).NoGc_half_sim.massElevator = Me_half;
            Systems(i).NoGc_half_sim.bCompensator = 0;
        Systems(i).NoGc_full_sim = inSimulink;
            Systems(i).NoGc_full_sim.shortName = ['NoGc-full-Sim'];
            Systems(i).NoGc_full_sim.longName = ['Closed-Loop Step:No Gc Full Load(SIMULINK)'];
            Systems(i).NoGc_full_sim.massElevator = Me_full;
            Systems(i).NoGc_full_sim.bCompensator = 0;

        if(200 == SimulationLevel)
            %Models(i).LagGc_empty_mat = f_Simulate(Systems(i).LagGc_empty_mat);
            Models(i).LagGc_half_mat = f_Simulate(Systems(i).LagGc_half_mat);
            Models(i).LagGc_full_mat = f_Simulate(Systems(i).LagGc_full_mat);
            %Models(i).NoGc_half_mat = f_Simulate(Systems(i).NoGc_full_mat);
            %Models(i).NoGc_full_mat = f_Simulate(Systems(i).NoGc_full_mat);

            %Models(i).LagGc_empty_sim = f_Simulate(Systems(i).LagGc_empty_sim);
            Models(i).LagGc_half_sim = f_Simulate(Systems(i).LagGc_half_sim);
            Models(i).LagGc_full_sim = f_Simulate(Systems(i).LagGc_full_sim);
            %Models(i).NoGc_half_sim = f_Simulate(Systems(i).NoGc_half_sim);
            %Models(i).NoGc_full_sim = f_Simulate(Systems(i).NoGc_full_sim);

            %Analyze = Models(1).LagGc_full_sim;
            %Models(2).LagGc_full_sim.Out.txtAnalysisSummary
        %Post analysis
            
            SystemsInModel = fields(Models(1));
            for j = 1:numel(Models)
            end
            
            figure();
            hold on
                plot(Models(i).LagGc_half_mat.Out.tout, Models(i).LagGc_half_mat.Out.Xx,'g');
                plot(Models(i).LagGc_full_mat.Out.tout, Models(i).LagGc_full_mat.Out.Xx,'r');
                plot(Models(i).LagGc_half_sim.Out.tout, Models(i).LagGc_half_sim.Out.Xx,'m');
                plot(Models(i).LagGc_full_sim.Out.tout, Models(i).LagGc_full_sim.Out.Xx,'b');
            hold off
                set(figure(gcf), 'WindowStyle', 'docked')
                title(['Closed-Loop Step:LagGc[',f_nameLagGc(Models(i).LagGc_full_sim.System.LagDesign),'] Position of Elevator']);
                legend(Models(i).LagGc_half_mat.System.shortName,Models(i).LagGc_full_mat.System.shortName,...
                    Models(i).LagGc_half_sim.System.shortName,Models(i).LagGc_full_sim.System.shortName);
        %% Plot Entire Model Xx with NoGc
%         figure();
%             hold on
%                 plot(Models(i).LagGc_half_mat.Out.tout, Models(i).LagGc_half_mat.Out.Xx,'g');
%                 plot(Models(i).LagGc_full_mat.Out.tout, Models(i).LagGc_full_mat.Out.Xx,'r');
%                 plot(Models(i).NoGc_full_mat.Out.tout, Models(i).NoGc_full_mat.Out.Xx,'.r');
%                 plot(Models(i).LagGc_half_sim.Out.tout, Models(i).LagGc_half_sim.Out.Xx,'m');
%                 plot(Models(i).LagGc_full_sim.Out.tout, Models(i).LagGc_full_sim.Out.Xx,'b');
%                 plot(Models(i).NoGc_full_sim.Out.tout, Models(i).NoGc_full_sim.Out.Xx,'.b');
%             hold off
%                 set(figure(gcf), 'WindowStyle', 'docked')
%                 title(['Closed-Loop Step:LagGc[',f_nameLagGc(Models(i).LagGc_full_sim.System.LagDesign),'] Position of Elevator']);
%                 legend(Models(i).LagGc_half_mat.System.shortName,Models(i).LagGc_full_mat.System.shortName,Models(i).NoGc_full_mat.System.shortName,...
%                     Models(i).LagGc_half_sim.System.shortName,Models(i).LagGc_full_sim.System.shortName,Models(i).NoGc_full_sim.System.shortName);

        end
    end
end
%% Iterative approach
if(SimulationLevel<0)
    SimulationLevel = abs(SimulationLevel);
    
    % Change these to iterate. const=-1 means it is the single varying variable.
    Iterate.k_const = -1;%LagDesigns(2).k;
    Iterate.Pc_const= StartingLagDesign.Pc;%0.0072;
    Iterate.a_const = StartingLagDesign.a;%LagDesigns(3).a;
    
    Iterate.k =     [250:25:375];%;linspace(250,375,SimulationLevel);
    %Iterate.k =     linspace(Iterate.k_const*0.8,Iterate.k_const*1.2,SimulationLevel);
    Iterate.Pc =    linspace(0.04,0.06,length(Iterate.k));
    Iterate.a =     linspace(1,5,SimulationLevel);
    %     Iterate.k =     linspace(100,300,SimulationLevel);
    %     Iterate.Pc =    linspace(0.001,1.5,SimulationLevel);
    %     Iterate.a =     linspace(1,25,SimulationLevel);

    inMatlab.bCompensator = 1; 
        inMatlab.bStepNotRamp = 1;
        inMatlab.bMatlabNotSimulink =1;   inMatlab.bGravity = 0;  inMatlab.bLimiters = 0;
        inMatlab.massElevator = Me_full; 
    inSimulink.bCompensator = 1; 
        inSimulink.bStepNotRamp = 1;
        inSimulink.bMatlabNotSimulink =0;   inSimulink.bGravity = 1;  inSimulink.bLimiters = 1;
        inSimulink.massElevator = Me_full;
    
    for iDesign3 = 1:1
            Iterate.Pc_const = Iterate.Pc_const;%Iterate.Pc(iDesign2);  
        for iDesign2 = 1:1
            Iterate.a_const = Iterate.a_const;%Iterate.a_const;%Iterate.a(iDesign3);
            
            if(Iterate.k_const==-1)
                Iterate.Constants = sprintf('(%g=a,%g=Pc)',Iterate.a_const,Iterate.Pc_const);
            elseif(Iterate.Pc_const==-1)
                Iterate.Constants = sprintf('(%g=a,%g=k)',Iterate.a_const,Iterate.k_const);
            elseif(Iterate.a_const==-1)
                Iterate.Constants = sprintf('(%g=Pc,%g=k)',Iterate.Pc_const,Iterate.k_const);
            end
            inMatlab.saveDir = ['EE315_Part3Results\IterateMatlab\',Iterate.Constants]; %Iterate.Constants;
            inSimulink.saveDir = ['EE315_Part3Results\IterateSimulink\',Iterate.Constants]; %Iterate.Constants;
            for iDesign = 1:SimulationLevel
                LagDesigns(iDesign).k = Iterate.k_const;
                LagDesigns(iDesign).Pc = Iterate.Pc_const;
                LagDesigns(iDesign).a = Iterate.a_const; %Iterate.a(iDesign);
                if(Iterate.k_const==-1)
                    LagDesigns(iDesign).k = Iterate.k(iDesign);
                elseif(Iterate.Pc_const==-1)
                    LagDesigns(iDesign).Pc = Iterate.Pc(iDesign);
                elseif(Iterate.a_const==-1)
                    LagDesigns(iDesign).a = Iterate.a(iDesign);
                end
            end
            for iDesign = 1:SimulationLevel
                if(Iterate.k_const==-1)
                    Iterate.Variable = sprintf('(%g=k)',Iterate.k(iDesign));
                elseif(Iterate.Pc_const==-1)
                    Iterate.Variable = sprintf('(%g=Pc)',Iterate.Pc(iDesign));
                elseif(Iterate.a_const==-1)
                    Iterate.Variable = sprintf('(%g=a)',Iterate.a(iDesign));
                end
                LagDesign = LagDesigns(iDesign);
                
                inMatlab.shortName = [Iterate.Variable,'LagGcNoGrav',num2str(iDesign),'-full-Mat'];
                    inMatlab.LagDesign = LagDesign; inMatlab.k = LagDesign.k;
                    inMatlab.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(inMatlab.LagDesign),'] No Gravity,Full Load (MATLAB)'];
                inSimulink.shortName = [Iterate.Variable,'LagGc',num2str(iDesign),'-full-Sim'];
                    inSimulink.LagDesign = LagDesign; inSimulink.k = LagDesign.k;
                    inSimulink.longName = ['Closed-Loop Step:LagGc[',f_nameLagGc(inSimulink.LagDesign),'] Full Load (SIMULINK)'];
                    
                inSimulink.Gc = f_Gc(inSimulink.LagDesign.Pc,inSimulink.LagDesign.a);
                    SimulinkGc.inMatlab(iDesign2*SimulationLevel+iDesign) = f_Simulate(inSimulink);
                inMatlab.Gc = f_Gc(inMatlab.LagDesign.Pc,inMatlab.LagDesign.a);
                    MatlabGc.inMatlab(iDesign2*SimulationLevel+iDesign) = f_Simulate(inMatlab);
            end
        end % Pc
    end % a
end



% E2out = f_Simulate(E2);
% E2_a = repmat(E2out,numel(Iterate.a),1);
% for i=1:1
%     sys = E2;
%     sys.a =Iterate.a(i);
%     %sys.k =
%     %sys.Pc
%     sys.shortName = [E2.shortName,'_', f_nameLagGc(sys)];
%     sys.Gc = f_Gc(sys.Pc,sys.a);
%     E2_a(i) = f_Simulate(sys);
% end

%clear Time Gc_num Gc_den K Kmd Height InputRamp_slope InputRamp_ivel SimulinkInputSelector UN Nm Ng Tg Deq J kt kb La Ra maxVoltageOfMotorDrive maxTorqueOfMotor maxTorqueOfMotor maxPosition r

%save([System.saveDir,'\1Workspace.mat'])




