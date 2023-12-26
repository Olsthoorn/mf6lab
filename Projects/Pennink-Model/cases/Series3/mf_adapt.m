%% Pennink1 - density flow in Pennnink's (1915) sand box model
% Experiments series 3
%
%see http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en
%
% The computation takes about 10 minutes. The computation time may be
% reduced to less than one minute at the cost of resolution. To do that,
% you may increase % the values of dx and dz from 1 cm to 2 cm below.
%
% in this experiment, Pennink studies simultaneous cFresh and saltwater
% flow in a cross section fed by rain and discharging to a canal at the
% westsided of the section.
% He first establiches a saltwater zone by letting milk enter slowly at the
% bottom of the model, regulated by a fixed-head reservoir.
% After the milk volum as reached 13 cm height in the model, he starts the
% precipitation 4.2 L/h over a width of 45 cm.
% This causes the interface to reach equilibrium.
% During this phase, milk is injected continuously at a rate of 2.14
% ml/min, the rate that fills the model to 13 cm in 5 hours time, as was
% the time Pennink took to let the milk enter the model. This rate seems
% also ok in the subsequent simulation. Pennink mentions that the milk was
% added drop by drop. THerefore, no fixed saltwater head was applied, but a
% well was used as the boundary condition for the milk entry point instead.
% After a day, Pennink stops adding milk and then after one more day, he
% takes the last photo of this test series.
% As it takes about two hours to establish virtual equilibrium of the
% saltwater interface, I shortened the simulation to 2 hours for the time
% to reach equilbrium with the milk injection and then continu during a
% full day without the milk injecition to simulate the washing out of the
% milk from the model. This should allow for calibration of the
% dispersivity in Penninks original model. This calibration was not really
% attempted as the situation after a full day without milk addition was
% similar enough the that of the last photo in Pennink's series 3 test.
% A somewhat better result may be obtained by raising the dispersivity a
% bit.
% Computation time about 30 min. To reduce the compuation time you may use
% a coarser grid, say 2 cm instead of 1 cm.
%
% As alway, the models can be readily run with mt3dms as with seawat by
% switching the appropriate packages in the NAM worksheet. An alternative
% to running mt3dms is to run SEAWAT without density active. This is done
% by setting -1 for the on switch for the VDF package in the NAM worksheet
% instead of 1.
%
% TO 090312 100523 100721

clear variables; close all;

basename='PenninkSeries3';
save name basename

GREP = 'STRESS PERIOD';

fprintf('Tests Pennink (1915), series 3\n');

xSand =[
   -2.2878
   67.6080
   67.7985
   43.8015
   30.0890
   23.2327
   17.7096
   11.4247
   10.4725
    8.5679
   -2.2878
   ];

zSand =[
   -2.3540
   -2.5445
   60.3133
   60.1229
   59.1705
   57.8371
   55.9323
   52.6942
   52.6942
   50.4085
   49.2656
   ];

xCanL =[
    -2.2878
    3.2353
    5.7112
    7.4252
    7.8061
   10.4725
   -2.2878
   ];

zCanL =[
   39.7417
   38.7893
   41.0750
   45.2656
   50.7894
   65.8372
   65.6467
   ];

xCanR=[
   68.3698
   66.2748
   61.1326
   61.3231
   61.8944
   63.0372
   65.5130
   69.1316
   68.7507
   ];

zCanR =[
   63.3610
   66.2182
   65.2658
   51.7418
   46.9799
   42.4084
   41.2655
   41.4560
   64.3134
];


%% The grid, the box is 96x96 cm and 1.8 cm thick
cFresh = 0;    % Relative minimum concentration
cSalt  = 1;    % Relative maximum concentration
k          = 86500/(24*60); % cm/min calibrated
peff       = 0.38;       % calibrated
ss         = 1e-4;   % storage coefficient

MW     = 65;         % [cm] Width of model. Pennink p6
MH     = 65;         % [cm] Height of mode.
D      =  2;         % [cm] thickness of model
Fringe = 2.0;        % [cm] Capillary fringe thicknesssee page 34
zIface = 13;         % [cm] inial interface elevation
hCanL  = 45;         % [cm] Left canal water elevation
hCanR  = 46;         % [cm] Right canal water elevation

milkInjPnt = [63.5 1 0.1]; % point of milk supply

%% zone numbers (arbitrary as long as not interfering with IBOUND)
iSand = 5;
iMilk = 6;
iCanL = 7;
iCanR = 8;

%% Grid (always 3D)
dx=1.0;              % [cm] grid cell width
dy=  D;              % [cm] grid cell length = depth of model
dz=1.0;              % [cm] grid cell size vertical

xGr=[0:dx:MW MW];    % [cm] grid, making sure that MW is always included irrespective of dx !
yGr=[0 D];           % [cm] grid
zGr=[MH:-dz:0 0];    % [cm] grid, making sure that 0 is always include irrespective of dz !

gr = gridObj(xGr,yGr,zGr);

idxMilk = xyzindex(milkInjPnt,gr);

%% Model arrays

IBOUND=gr.const(0);

% implant zone numbers
IBOUND(inpolyz(gr.XM,gr.ZM,xSand,zSand))=iSand;
IBOUND(inpolyz(gr.XM,gr.ZM,xCanL,zCanL))=iCanL;
IBOUND(inpolyz(gr.XM,gr.ZM,xCanR,zCanR))=iCanR;

IBOUND(IBOUND==iCanL & gr.ZM>hCanL) = 0;  % dry make inactive
IBOUND(IBOUND==iCanR & gr.ZM>hCanR) = 0;  % dry make inactive

IBOUND(idxMilk)=iMilk;  % milk injection point

HK    = gr.const(k);
VK    = gr.const(k);
HK(HK>0 & gr.ZM>hCanL+Fringe)=k/10; % HK=0 above full capillary zone

SS    = gr.const(ss);  % storage/transient damps density change shocks

PEFF  = gr.const(peff);
PEFF(gr.ZM>hCanL+Fringe)=peff/3;  %unsaturated

STRTHD = gr.const(mean(hCanL,hCanR));

STCONC               = gr.const(cFresh);
STCONC(gr.ZM<zIface) = cSalt; 

ICBUND = IBOUND;

%% Get period data to set CHD boundaries for MODFLOW and PNTSRC for MTRDMS/SSM

zoneVals1 = {iCanL,'hL' 'hL'};
zoneVals2 = {iMilk,'QMilk'};
    
NPER= getPeriods(basename);
[CHD,PNSTRC              ] = bcnZone(basename,'CHD',IBOUND, zoneVals1,{cFresh});
%[WEL,PNSTRC(end+(1:NPER))] = bcnZone(basename,'WEL',IBOUND, zoneVals2,{cSalt});

%% RCH Recharge see below N=4.2 L/h according to Pennink

W=45; % [cm] width of rain added to top of model
M=37; % [cm] center of rain added to top of model (mid between 2 screws see photo)
N=4.2*1e3/60/W/dy;  % recharge in cm/min: 4.2 L/h / 60 min/h / W / dy

RECH =(gr.Xm>=M-0.5*W & gr.Xm<=M+0.5*W) * N;
RECH = bsxfun(@times,RECH,YS(1:NPER));

save underneath xSand zSand xCanL zCanL % needed to set PSIMASK
