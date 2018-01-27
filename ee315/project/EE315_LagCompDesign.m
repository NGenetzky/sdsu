clear all
close all
clc

f_moveFigure = @(H) set(figure(H), 'WindowStyle', 'docked');
Type        = 0;    % Type = [0,1,2,3] = [all,a,K,Pc]
ClassDesign = 1;    % ClassDesign = [1,~1]
%% setup variables
Jm = 11;                            % Kg*m^2
Dm = 0.117;                         % N*m*s*(rad^-1)
Ra = 54.8 *10^-3;                   % Ohms
La = 0.92 *10^-3;                   % H
UN = 470;                           % Vdc
kt = 3285/675;                      % V*(rad^-1)*s  
kb = (UN-Ra*675)/(825*(2*pi/60));   % N*(m^-1) %kb=(UN-Ra*IN)/(n(2pi/60)) % n=825

Jg        = 250;      % Kg*m^2
MeEMPTY   = 7500;     % Kg
r         = 1;        % m
Dg        = 75.8;     % N*m*s*(rad^-1)
De        = 151.6;    % N*s*(m^-1) 
Nm        = 60;       % Gear Teeth 
Ng        = 1440;     % Gear Teeth

MassPeople    = 62 * 10;     %Average mass of Human * Average number of people
Me            = MeEMPTY + MassPeople;
Deq           = Dm + 4*Dg*(Nm/Ng)^2+2*De*r^2*(Nm/Ng)^2;
J             = Jm + 4*Jg*(Nm/Ng)^2+Me*r^2*(Nm/Ng)^2;
Ag            = 9.82;
Tg            = (Nm/Ng)*r*Ag*Me;

Amp       = 1;    % r(s) = 1/s
Fudge     = 9;    % 'Fudge' factor [rad/s]
Height = 50;
w_log   = logspace(-2,3,1000);
%% setup initial TF  
Gx.b0 = kt*r*(Nm/Ng)*1/(La*J);
Gx.a0 = 0;
Gx.a1 = (Deq*Ra+kt*kb)/(La*J);
Gx.a2 = (Ra*J+Deq*La)/(La*J);
Gx.a3 = 1;
    Gx.num  = Gx.b0;
    Gx.den  = [Gx.a3 Gx.a2 Gx.a1 Gx.a0];
    Gx.TF   = tf(Gx.num,Gx.den);
        figure
        bode(Gx.TF,w_log),hold on,grid on,f_moveFigure(gcf);
        title('Bode of Transfer Function G');
        G = Gx.TF
Height    = 50;


    
    
opt = stepDataOptions('StepAmplitude',Height);
%% 
if(ClassDesign)
    Desired.OS      = 12;
    Desired.e_ss    = 0.1;
        fprintf('Initial constraints are a %%O.S. of %0.4f %% and Steady State Error of %0.4f\n', Desired.OS,Desired.e_ss)
    Desired.Zeta    = -log(Desired.OS/100)/sqrt(pi^2+(log(Desired.OS/100))^2);
    Desired.PM      = atand(2*Desired.Zeta/sqrt((-2*Desired.Zeta^2+sqrt(1+4*Desired.Zeta^4))));
        fprintf('Desired PM for %%OS is %6.4f degrees\n',Desired.PM)
    Desired.Kp      = Amp/Desired.e_ss;
        fprintf('Required Kp for Desired error is %6.4f\n',Desired.Kp)
    Gc.num          = Gx.num;
    Gc.den          = [Gx.a3 Gx.a2 Gx.a1];
    Gc.Kp           = polyval(Gc.num,0)/polyval(Gc.den,0);
        fprintf('Actual Kp for Initial TF is %6.4f\n',Gc.Kp)
    Gc.K            = Desired.Kp/Gc.Kp;
        fprintf('Required Gain K is %6.4f\n',Gc.K)
    Gc.KG           = Gc.K*Gx.TF;
        figure
        bode(Gc.KG,w_log),hold on,grid on,f_moveFigure(gcf);
        title('Bode of Transfer Function KG');
        KG = Gc.KG
    Desired.Phase           = -180+Desired.PM+Fudge;
        fprintf('The phase at 0 [dB] is %6.4f degrees\n',Desired.Phase)
    [Gc.KG_mag,Gc.KG_phase] = bode(Gc.KG,w_log);
    Gc.KG_mag               = Gc.KG_mag(1,:,:);
    [ans,Gc.WpmIND]         = min(abs(Gc.KG_phase(1,1,:)-(Desired.Phase) ));
    Gc.Wpm                  = w_log(Gc.WpmIND);
        fprintf('Which means the frequency of PM is %6.4f [rad/s]\n',Gc.Wpm)
    Gc.Pc                   = Gc.Wpm/100;
        fprintf('Therefore the pole is located at %6.4f, two decades before Wpm\n',Gc.Pc)
    Gc.Gain                 = db(Gc.KG_mag(Gc.WpmIND));
        fprintf('Current gain is %6.4f [dB]\n',Gc.Gain)
    
    a=[1:0.2:3];
    pc = 1;
    for k=1:length(a)
        legendInfo{k} = ['a = ' num2str(a(k))];
        zc=a(k)*pc;
        n=[1/zc 1];
        d=[1/pc 1];
        Gc.test=tf(n,d);
        figure(3)
        bode(Gc.test,w_log),hold on, grid on
        title('Gain of Normalized TF with varying ''a'''),legend(legendInfo), f_moveFigure(gcf)
    end
    Gc.a    = 1.4;
    fprintf('What value of ''a'' gets a gain of -%6.5f [dB] at %6.5f [rad/s]: %6.4f\n',Gc.Gain,Gc.Wpm,Gc.a)
    fprintf('Value of ''a'' chosen is %6.4f\n', Gc.a)
    Gc.Zc   = Gc.a*Gc.Pc;
    fprintf('Therfore a zero is located at %6.4f', Gc.Zc)
    Gc.num  = [1/Gc.Zc 1];
    Gc.den  = [1/Gc.Pc 1];
    Gc.TF   = tf(Gc.num,Gc.den);
    Gc.KG   = Gc.K*Gx.TF;
    KGcG    = series(Gc.TF,Gc.KG);
    T       = feedback(KGcG,1);
        
        figure
        bode(Gc.TF,w_log),hold on,grid on,f_moveFigure(gcf);
        title('Bode of Transfer Function Gc');
        Gc = Gc.TF
        figure
        bode(KGcG,w_log),hold on,grid on,f_moveFigure(gcf);
        title('Bode of Transfer Function KGcG');
        KGcG = KGcG
        figure
        step(T,opt),hold on,grid on,f_moveFigure(gcf);
        title('Plot of Closed-Loop Tranfer Function T(s) with a step input applied');
        T = T
else
switch(Type)
    case 0  % 3-D iterations
        Gc.K    = linspace(100,500,5);  %linspace(100,1000,5); %linspace(100,7000,5); %[200:50:500];
        Gc.Pc   = linspace(0.02,0.1,5); %linspace(0.01,1,5);   %linspace(0.001,1,5);  %[0.02:0.02:0.1];
        Gc.a    = linspace(3,8,5);      %linspace(2,10,8);     %linspace(2,10,5);     %[3:1:8];
        for i = 1:length(Gc.K)
           for j = 1:length(Gc.Pc)
               titleInfo{i,j} = ['K = ' num2str(Gc.K(i)), '  Pc = ' num2str(Gc.Pc(j))];
               for k = 1:length(Gc.a)
                    legendInfo{k}   = ['a = ' num2str(Gc.a(k))];
                    Gc.Zc(j,k)    = Gc.a(k)*Gc.Pc(j);
                    Gc.num          = [1/Gc.Zc(j,k) 1];
                    Gc.den          = [1/Gc.Pc(j) 1];
                    Gc.TF           = tf(Gc.num,Gc.den);
                    Gc.KG           = Gc.K(i)*Gx.TF;
                    KGcG            = series(Gc.TF,Gc.KG);
                    T               = feedback(KGcG,1);
                        figure(10*i+j),step(T,opt),hold on,grid on;
                        title(titleInfo(i,j)),legend(legendInfo),f_moveFigure(gcf);
               end        
            end
        end
    case 1  % 1-D varying a
        Gc.K        = Desired.Kp/Gc.Kp;      
        Gc.Pc       = Gc.Wpm/100;
        Gc.a        = linspace(2,11,10);
        titleInfo   = ['K = ' num2str(Gc.K), '  Pc = ' num2str(Gc.Pc)];
        for k = 1:length(Gc.a)
            legendInfo{k}   = ['a = ' num2str(Gc.a(k))];
            Gc.Zc(j,k)    = Gc.a(k)*Gc.Pc(j);
            Gc.num          = [1/Gc.Zc(j,k) 1];
            Gc.den          = [1/Gc.Pc(j) 1];
            Gc.TF           = tf(Gc.num,Gc.den);
            Gc.KG           = Gc.K(i)*Gx.TF;
            KGcG            = series(Gc.TF,Gc.KG);
            T               = feedback(KGcG,1);
                figure(1),step(T,opt),hold on,grid on;
                title(titleInfo),legend(legendInfo),f_moveFigure(gcf);
        end        
    case 2  % 1-D varying K
        Gc.K        = [225:25:375];     %linspace(100,1000,10);      
        Gc.Pc       = 0.05;             %Gc.Wpm/100;
        Gc.a        = 4.6;
        titleInfo   = ['Pc = ' num2str(Gc.Pc), '  a = ' num2str(Gc.a)];
        for k = 1:length(Gc.K)
            legendInfo{k}   = ['K = ' num2str(Gc.K(k))];
            Gc.Zc           = Gc.a*Gc.Pc;
            Gc.num          = [1/Gc.Zc 1];
            Gc.den          = [1/Gc.Pc 1];
            Gc.TF           = tf(Gc.num,Gc.den);
            Gc.KG           = Gc.K(k)*Gx.TF;
            KGcG            = series(Gc.TF,Gc.KG);
            T               = feedback(KGcG,1);
                figure(1),step(T,opt),hold on,grid on;
                title(titleInfo),legend(legendInfo),f_moveFigure(gcf);
        end 
    case 3  % 1-D varying Pc
        Gc.K        = Desired.Kp/Gc.Kp;     
        Gc.Pc       = linspace(0.001,1,10);
        Gc.a        = 4;
        titleInfo   = ['K = ' num2str(Gc.K), '  a = ' num2str(Gc.a)];
        for k = 1:length(Gc.Pc)
            legendInfo{k}   = ['Pc = ' num2str(Gc.Pc(k))];
            Gc.Zc           = Gc.a*Gc.Pc(k);
            Gc.num          = [1/Gc.Zc 1];
            Gc.den          = [1/Gc.Pc(k) 1];
            Gc.TF           = tf(Gc.num,Gc.den);
            Gc.KG           = Gc.K*Gx.TF;
            KGcG            = series(Gc.TF,Gc.KG);
            T               = feedback(KGcG,1);
                figure(1),step(T,opt),hold on,grid on;
                title(titleInfo),legend(legendInfo),f_moveFigure(gcf);
        end       
end
end
    %%
    