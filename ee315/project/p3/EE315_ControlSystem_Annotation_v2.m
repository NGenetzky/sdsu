function [ Annote ] = EE315_ControlSystem_Annotation_v2(System,Time,Type)
%% System variables
% System.Overshoot
% System.PeakPosition
% System.PeakTimePosition
% System.RiseTimePosition
% System.SteadyPosition
% System.SettlingTimePosition
System.PeakOffsetScale = 0.03;
%% Units for variables
Figures.Units.distance      = '[m]';
Figures.Units.TimePosition  = '[s]';
Figures.Units.percent       = '%';
Figures.Units.velocity      = '[m/s]';
Figures.Units.torque        = '[N*m]';
Figures.Units.current       = '[A]';
%% annotation values
Annote = 1;
switch(Type)
    case 1 % position
        System.loc.Xx.PeakOffset   = System.PeakTime+System.PeakTime*System.PeakOffsetScale;
        System.loc.Xx.Peak         = System.PeakTime;
        System.loc.Xx.Rise         = System.RiseTime; 
        System.loc.Xx.Settling     = System.SettlingTime;
        System.loc.Xx.Overshoot     = 0.5;

        System.loc.Xy.PeakOffset   = System.PeakPosition-System.PeakPosition*System.PeakOffsetScale;
        System.loc.Xy.Peak         = System.PeakPosition;
        System.loc.Xy.Rise         = System.RisePosition;
        System.loc.Xy.Steady       = System.SteadyPosition;
        System.loc.Xy.Settling     = System.SettlingPosition; 
        System.loc.Xy.Overshoot    = (System.SteadyPosition+System.PeakPosition)/2;
        System.loc.Xy.UpperTol     = System.SteadyPosition+System.SteadyPosition*System.ToleranceForTsTr;
        System.loc.Xy.LowerTol     = System.SteadyPosition-System.SteadyPosition*System.ToleranceForTsTr;
        % annotation plots/texts
        hold on
        line(System.loc.Xx.Peak,                                System.loc.Xy.Peak,                             'Color',[0 204 0]/255,      'LineStyle','.','MarkerSize',15); % point for PeakPosition
        line([System.loc.Xx.Peak,System.loc.Xx.PeakOffset],     [System.loc.Xy.Peak,System.loc.Xy.PeakOffset],  'Color',[0 204 0]/255);                         % line for PeakPosition
        line([0,System.loc.Xx.Peak],                            [System.loc.Xy.Peak,System.loc.Xy.Peak],        'Color',[0 204 0]/255,      'LineStyle','--');  % horizontal line for PeakPosition
        line([System.loc.Xx.Peak,System.loc.Xx.Peak],           [0,System.loc.Xy.Peak],                         'Color',[0 204 0]/255,      'LineStyle','--');  % vertical line for PeakPosition
        line([System.loc.Xx.Rise,System.loc.Xx.Rise],           [0,System.loc.Xy.Rise],                         'Color',[255 153 51]/255,   'LineStyle','--');  % vertical line for RiseTimePosition
        line([0,Time.Stop],                                     [System.loc.Xy.Steady,System.loc.Xy.Steady],    'Color',[255 0 0]/255,      'LineStyle','--');  % horizontal line for SteadyPosition
        line([0,Time.Stop],                                     [System.loc.Xy.UpperTol,System.loc.Xy.UpperTol],'Color',[204 204 0]/255,    'LineStyle','--');  % horizontal line for SteadyPosition UpperTol
        line([0,Time.Stop],                                     [System.loc.Xy.LowerTol,System.loc.Xy.LowerTol],'Color',[204 204 0]/255,    'LineStyle','--');  % horizontal line for SteadyPosition LowerTol
        line([System.loc.Xx.Settling,System.loc.Xx.Settling],   [0,System.loc.Xy.Settling],                     'Color',[255 153 51]/255,   'LineStyle','--');  % vertical line for SettlingTimePosition
        line([System.loc.Xx.Overshoot,System.loc.Xx.Overshoot], [System.loc.Xy.Steady,System.loc.Xy.Peak],      'Color',[0 0 0]/255,        'LineStyle',':');   % vertical line for Overshoot
        %------------------------------------------------------------------------------------------------------------------------------- 
        text(System.loc.Xx.PeakOffset,System.loc.Xy.PeakOffset, ['Peak Position = ',                    num2str(System.PeakPosition),           ' ',Figures.Units.distance],        'VerticalAlignment','top',      'HorizontalAlignment','left','color',[0 204 0]/255,     'FontSize',8);  % label for PeakPosition
        text(0.001,System.loc.Xy.Peak,                          ['Peak Position = ',                    num2str(System.PeakPosition),           ' ',Figures.Units.distance],        'VerticalAlignment','top',      'HorizontalAlignment','left','color',[0 204 0]/255,     'FontSize',8);  % label for horizontal line for PeakPosition
        text(System.loc.Xx.Peak,2,                              ['  Time to Peak Position = ',          num2str(System.PeakTime),               ' ',Figures.Units.TimePosition],    'VerticalAlignment','bottom',   'HorizontalAlignment','left','color',[0 204 0]/255,     'FontSize',8);  % label for vertical line for PeakPosition
        text(0.001,System.loc.Xy.Steady,                        ['Steady State Position = ',            num2str(System.SteadyPosition),         ' ',Figures.Units.distance],        'VerticalAlignment','top',      'HorizontalAlignment','left','color',[255 0 0]/255,     'FontSize',8);  % label for horizontal line for SteadyPosition
        text(0.001,System.loc.Xy.LowerTol-0.5,                  ['Steady State Position Range: \pm ',   num2str(System.ToleranceForTsTr*100),                ' ',Figures.Units.percent],         'VerticalAlignment','top',      'HorizontalAlignment','left','color',[204 204 0]/255,   'FontSize',8);  % label for horizontal lines for Tolerance Range
        text(System.loc.Xx.Rise,2,                              ['  RiseTime = ',                       num2str(System.RiseTime),               ' ',Figures.Units.TimePosition],    'VerticalAlignment','bottom',   'HorizontalAlignment','left','color',[255 153 51]/255,  'FontSize',8);  % label for vertical line for RiseTimePosition
        text(System.loc.Xx.Settling,1,                          ['  Settling Time = ',                  num2str(System.SettlingTime),           ' ',Figures.Units.TimePosition],    'VerticalAlignment','bottom',   'HorizontalAlignment','left','color',[255 153 51]/255,  'FontSize',8);  % label for vertical line for SettlingTimePosition
        text(System.loc.Xx.Overshoot,System.loc.Xy.Overshoot,   ['  % Overshoot = ',                    num2str(System.Overshoot),              ' ',Figures.Units.percent],         'VerticalAlignment','cap',      'HorizontalAlignment','left','color',[0 0 0]/255,       'FontSize',8);  % label for vertical line for Overshoot
    case 2 % torque
        Figures.Annotate.TmLoc.xPeakOffset  = System.PeakTimeMotorTorque+System.PeakTimeMotorTorque*1;
        Figures.Annotate.TmLoc.xPeak        = System.PeakTimeMotorTorque;
        Figures.Annotate.TmLoc.xSettling    = System.SteadyTimeMotorTorque;
        Figures.Annotate.TmLoc.yPeakOffset  = System.PeakMotorTorque-System.PeakMotorTorque*0.03;
        Figures.Annotate.TmLoc.yPeak        = System.PeakMotorTorque;
        Figures.Annotate.TmLoc.ySteady      = System.SteadyMotorTorque;
        Figures.Annotate.TmLoc.yOvershoot   = (System.SteadyMotorTorque+System.PeakMotorTorque)/2;
        Figures.Annotate.TmLoc.yUpperTol    = System.SteadyMotorTorque+System.SteadyMotorTorque*System.ToleranceForTsTr;
        Figures.Annotate.TmLoc.yLowerTol    = System.SteadyMotorTorque-System.SteadyMotorTorque*System.ToleranceForTsTr;
        hold on 
        line(Figures.Annotate.TmLoc.xPeak,                                       Figures.Annotate.TmLoc.yPeak,                                      'Color',[0 204 0]/255,  'LineStyle','.','MarkerSize',15);   % point for PeakTorque
        line([Figures.Annotate.TmLoc.xPeakOffset,Figures.Annotate.TmLoc.xPeak], [Figures.Annotate.TmLoc.yPeakOffset,Figures.Annotate.TmLoc.yPeak],  'Color',[0 204 0]/255);                                     % line for PeakTorque
        line([0,Figures.Annotate.TmLoc.xPeak],                                  [Figures.Annotate.TmLoc.yPeak,Figures.Annotate.TmLoc.yPeak],        'Color',[0 204 0]/255,  'LineStyle','--');                  % horizontal line for PeakTorque
        line([Figures.Annotate.TmLoc.xPeak,Figures.Annotate.TmLoc.xPeak],       [0,Figures.Annotate.TmLoc.yPeak],                                   'Color',[0 204 0]/255,  'LineStyle','--');                  % vertical line for PeakTorque
        line([0,Time.Stop],                                                     [Figures.Annotate.TmLoc.ySteady,Figures.Annotate.TmLoc.ySteady],    'Color',[255 0 0]/255,  'LineStyle','--');                  % horizontal line for SteadyTorque
        %------------------------------------------------------------------------------------------------------------------------------- 
        text(Figures.Annotate.TmLoc.xPeakOffset,    Figures.Annotate.TmLoc.yPeakOffset, ['Peak Torque = ',          num2str(System.PeakMotorTorque),   ' ',Figures.Units.torque],'VerticalAlignment','top',    'HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8); % label for PeakTorque
        text(0,                                     Figures.Annotate.TmLoc.ySteady,     ['Steady State Torque = ',  num2str(System.SteadyMotorTorque), ' ',Figures.Units.torque],'VerticalAlignment','bottom', 'HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8); % label for horizontal line for SteadyTorque
    case 3 % current
        Figures.Annotate.IaLoc.xPeakOffset  = System.PeakTimeMotorCurrent+System.PeakTimeMotorCurrent*1;
        Figures.Annotate.IaLoc.xPeak        = System.PeakTimeMotorCurrent;
        Figures.Annotate.IaLoc.xSettling    = System.SteadyTimeMotorCurrent;
        Figures.Annotate.IaLoc.yPeakOffset  = System.PeakMotorCurrent-System.PeakMotorCurrent*0.03;
        Figures.Annotate.IaLoc.yPeak        = System.PeakMotorCurrent;
        Figures.Annotate.IaLoc.ySteady      = System.SteadyMotorCurrent;
        Figures.Annotate.IaLoc.yOvershoot   = (System.SteadyMotorCurrent+System.PeakMotorCurrent)/2;
        Figures.Annotate.IaLoc.yUpperTol    = System.SteadyMotorCurrent+System.SteadyMotorCurrent*System.ToleranceForTsTr;
        Figures.Annotate.IaLoc.yLowerTol    = System.SteadyMotorCurrent-System.SteadyMotorCurrent*System.ToleranceForTsTr;
        hold on 
        line(Figures.Annotate.IaLoc.xPeak,                                       Figures.Annotate.IaLoc.yPeak,                                      'Color',[0 204 0]/255,  'LineStyle','.','MarkerSize',15);   % point for PeakCurrent
        line([Figures.Annotate.IaLoc.xPeakOffset,Figures.Annotate.IaLoc.xPeak], [Figures.Annotate.IaLoc.yPeakOffset,Figures.Annotate.IaLoc.yPeak],  'Color',[0 204 0]/255);                                     % line for PeakCurrent
        line([0,Figures.Annotate.IaLoc.xPeak],                                  [Figures.Annotate.IaLoc.yPeak,Figures.Annotate.IaLoc.yPeak],        'Color',[0 204 0]/255,  'LineStyle','--');                  % horizontal line for PeakCurrent
        line([Figures.Annotate.IaLoc.xPeak,Figures.Annotate.IaLoc.xPeak],       [0,Figures.Annotate.IaLoc.yPeak],                                   'Color',[0 204 0]/255,  'LineStyle','--');                  % vertical line for PeakCurrent
        line([0,Time.Stop],                                                     [Figures.Annotate.IaLoc.ySteady,Figures.Annotate.IaLoc.ySteady],    'Color',[255 0 0]/255,  'LineStyle','--');                  % horizontal line for SteadyCurrent
        %------------------------------------------------------------------------------------------------------------------------------- 
        text(Figures.Annotate.IaLoc.xPeakOffset,    Figures.Annotate.IaLoc.yPeakOffset, ['Peak Current = ',         num2str(System.PeakMotorCurrent),  ' ',Figures.Units.current],'VerticalAlignment','top',   'HorizontalAlignment','left','color',[0 204 0]/255,'FontSize',8); % label for PeakCurrent
        text(0,                                     Figures.Annotate.IaLoc.ySteady,     ['Steady State Current = ', num2str(System.SteadyMotorCurrent),' ',Figures.Units.current],'VerticalAlignment','bottom','HorizontalAlignment','left','color',[255 0 0]/255,'FontSize',8); % label for horizontal line for SteadyCurrent
end
end
