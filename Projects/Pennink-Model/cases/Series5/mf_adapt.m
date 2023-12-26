%% Pennink Series 5 - Marck 1905
% -- FIFTH SERIES OF EXPERIMENTS (MARCH 1905)
%
%see http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en
%
% in this experiment, Pennink studies simultaneous cFresh and saltwater
% flow from a canal in the right and one in the left corner of the sand-box model
% to a canal in the center. The following text is cited from his book, p63:
%
% Pennink writes:
% "To accelerate the motion and to observe the phenomena more quickly and
% more clearly, the slopes are now made steeper. For this purpose a canal
% has been reserved in the middle of hte sand-mass in which a lower level
% may be kept up.
% The discharge from that canal is brought about by means of a pipe
% connected with the back part of the apparatus, the cock of which is
% visible halfway on the left side."
%
% Pennink simulates upconing below the central canal. Interstingly he also
% studied upconing as influenced by a large lateral groundwater flow. This
% gives an interesting shape of the cone. It shows that it is possible for
% the water below the canal to remain fresh, while the salt water enters
% from the side.
%
% You may reduce the computation time by makeing the mesh coarser, say 2 cm
% instead of 1 cm.
%
% TO 090320 100508 100719

clear variables; close all;

basename='PenninkSeries5'; fprintf('Model %s\n',basename);

GREP = 'STRESS PERIOD';

%% Parameters
delta  = 0.03;   % Pennink used milk with density 1.028 - 1.03 kg/L, p48
cFresh = 0;      % Relative conc of cFresh
cSalt  = 1;      % Relative conc of saltwater (milk)
peff   = 0.38;   % [-] calibrated effective porosity
k   =86500/(24*60);       % m/min calibrated conductivity

MW    = 96;      % [cm] width  of model Pennink p63
MH    = 97;      % [cm] height of model Pennink p63
D     = 1.8;     % [cm] thickness of model

% All coordinates in cm relative to LL of sand mass in model

% countor around sand mass
xSand =[
   -3.6406
   97.8254
   97.8254
   46.8008
   -4.2237
];

zSand =[
   -0.7948
   -2.4891
   92.1082
   90.6963
   91.2611
];

% Line along top of sand mass
xTop =[   -5.5143    46.9644    97.7441];
zTop =[   90.6709    90.4615    91.7179];

% countour around right canal
xCanR =[
  100.5756
   86.6065
   90.3819
   92.0809
   94.9125
   96.9890
];

zCanR =[
   97.4417
   97.0229
   87.2506
   82.6437
   80.9685
   80.8987
];
   
% contour around left canal
xCanL =[
   -4.7592
    0.9040
    2.4142
    3.5468
    5.6233
   -4.1929
];

zCanL =[
   81.1081
   81.0383
   81.8061
   86.9016
   96.3947
   96.3947
];

% contour around center (mid) canal
xCanM =[
   41.3012
   43.7552
   44.5103
   45.8317
   48.2858
   49.4184
   50.1735
   51.1174
   51.4949
];

zCanM =[
   96.2551
   80.9685
   76.2219
   74.4769
   74.4769
   75.3145
   77.9670
   88.9957
   97.0229
];

zIFace=23; % initial elevation of interface

milkPoint = [91.5 1 1.5];  % coordates of point where milk was supplied

phiR = 87.5;               % [cm] head in right hand canal
phiS = 87.5-6.5;           % [cm] saline had as mentione by Pennink
Q    = 60e3/60;            % [cm/min, Pennink gives 60 L/h

% Elevation of water table including level in canals
xWT =[
   -3.5555
    3.3253
   44.4849
   50.2397
   90.5236
   98.6555
];

zWT =[
   87.4433
   87.3880
   77.3756
   77.3756
   87.2220
   87.3327
];

%% Model grid
dx = 2.0;              % [cm] grid horizontal cell size
dy =  D;               % [cm] thickness of model as specified by Pennink
dz = 2.0;              % [cm] grid vertical cell size

xGr = 0:dx:MW;         % xGr-grid coordinates
yGr = [0 dy];          % we take the model to be 1 m wide and adjust Pennink's discharge to it
zGr = (MH:-dz:0)';     % zGr line elevations

gr = gridObj(xGr,yGr,zGr);

idxMilk = xyzindex(milkPoint,gr);

ZT = XS( bsxfun(@times, gr.zm(:),interp1(xTop,zTop,gr.xm)));
WT = XS( bsxfun(@times, gr.zm(:),interp1(xWT ,zWT ,gr.xm)));

PhiS = WT(end)-6.5;  % [cm] salwater head according to Pennink

IBOUND = gr.const(0);

iSand = 1;
iMilk = 3;
iCanL = 4;
iCanM = 5;
iCanR = 6;

IBOUND(inpolyz(gr.XM,gr.ZM,xSand,zSand)) = iSand;
IBOUND(inpolyz(gr.XM,gr.ZM,xCanL,zCanL) & gr.ZM<WT) = iCanL;
IBOUND(inpolyz(gr.XM,gr.ZM,xCanM,zCanM) & gr.ZM<WT) = iCanM;
IBOUND(inpolyz(gr.XM,gr.ZM,xCanR,zCanR) & gr.ZM<WT) = iCanR;

IBOUND(idxMilk) = iMilk;

HK = gr.const(k); HK(gr.ZM>WT)=0.1*k;
VK = HK;

PEFF= gr.const(peff); PEFF(PEFF>WT) = peff/3;

STRTHD=WT; STRTHD(gr.ZM<zIFace)=PhiS;

ICBUND = IBOUND;

STCONC = gr.const(cFresh);      % Start with all cells no cSalt
STCONC(gr.ZM<zIFace) = cSalt;   % Make cells below interface salt

%% Head boundary conditions
% Pennink used head boundary conditions in the "canals" of his model. He
% changed these on several moments during the tests. These boundary
% condition are implemented in this simulation in using the CHD package.
% The values are in the accompanying workbook in sheet CHD
% As are the times that the photos were taken.

zoneVals = {iCanL 'hL' 'hL';
            iCanM 'hM' 'hM';
            iCanR 'hR' 'hR';
            iMilk 'hS' 'hS'};
        
zoneConc = {cFresh; cFresh; cFresh; cSalt};
        
[CHD,PNTSRC] = bcnZone(basename,'CHD',IBOUND,zoneVals,zoneConc);

save underneath xSand zSand xCanL zCanL xCanR zCanR xCanM zCanM % need for spi mask
