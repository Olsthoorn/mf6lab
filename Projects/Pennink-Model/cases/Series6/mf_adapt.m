%% Pennink Series 6 - Density flow simulation of Pennnink's (1915) sand-box model tests
% -- SIXTH SERIES OF EXPERIMENTS (MARCH 1905)
%
%see http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en
%
% in this experiment, Pennink studies simultaneous freshwater and saltwater
% flow from a canal in the right and one in the left corner of the sand-box model
% to a canal in the center. The following text is cited from his book, p87:
%
% "C. Experiments like sub B, but central-canal as a means of concentration
% substituted by a well-pipe with gauze strainer, placed in the canal axis.
%  Sixth series of experiments (March 1905)
% The water is discharged by means of the well-pipe, by siphoning, on the
% left side of the apparatus. Moreover, glass gauging-pipes are place in
% the sand in order to determine the pressure of the water just above the
% milk, and beside the gause strainer."
%
% The exact construction of the "well-pipe" is not discribed by Pennink, neither  
% nor its depth and screen length. In the simulatons done here, a fully screened
% well is assumed to a depth decuced from the photos in the book.
%
% We only model the situaton after Pennink pulled up his well screen to its
% final position of 32 cm above the original interface.
%
% Erroneous milk-head values may cause the milk to flow fast into or out of the model.
%
% Measuring directly from the photo was difficult and computing was not accurate
% enough due to vertical gradients in the sand above the interface.
% Therefore, the head in the milk reservoir was fine-tuned in several trial runs.
% The applied head values are in the workbook on sheet PER.

% TO 090320 100520

clear variables; close all;

basename='PenninkSeries6'; fprintf('Model %s\n',basename);

GREP='STRESS PERIOD';
BACKGROUND = false;

%% Parameters
Cfresh  = 0;             % Relative concentration of freshwater
cSalt    = 1;             % Relative concentration of saltwater (i.e. the milk)
peff    = 0.38;          % [-] calibrated effective porosity
k       = 86500/(24*60); % [cm/min] calibrated conductivity
GHBCOND = 1e3;           % [cm3/min] conductance of GHBs

%% All coordinates relative to the LL corner of the sand mass in the model
MH      = 96;        % [cm] Model height. Pennink specified this on p63
MW      = 97;        % [cm] Model width.  Pennink specified this on p63
D       = 1.8;       % [cm] Model thickness
hCanal0 = 91;        % [cm] measured from photo p88
zIFace  = 23;        % [cm]
zScrBot = zIFace+32; % [cm] 32 cm above original milk surface
zScrTop = MH;
%% Zone numbers

iUnsat= 1;
iSat  = 2;
iCanL = 3;
iWatL = 4;
iCanM = 5;
iWatM = 6;
iCanR = 7;
iWatR = 8;
iMilk = 9;
iWell = 10;

%% Zone contours

% Countour around the entire sand mass in the model
xSand = [   -0.0   97.0   97.0   48.5   -0.0 ];
zSand = [   -0.0   -0.0   92.1   90.7   92.1 ];

% Contour around the saturated sand body
xSat  = [   -0.0   97.0   97.0   48.5   -0.0 ];
zSat  = [   -0.0   -0.0   88.0   85.0   88.0 ];

% Contour around the unsaturated sand body
xUnsat  = [   -0.0   97.0   97.0   48.5   -0.0 ];
zUnsat  = [   97.0   97.0   88.0   85.0   88.0 ];

% Contour around left canal
xCanL = [   -0.0    0.9    1.5    2.0    2.5   -0.0 ];
zCanL = [   81.0   81.0   81.8   88.0   96.0   96.0 ];

% Water in the left canal
xWatL = [   -0.0    0.9    1.5    2.0   -0.0 ];
zWatL = [   81.0   81.0   81.8   88.0   88.0 ];

xCanR = MW-xCanL;
zCanR = zCanL;

% Water in the right-hand canal
xWatR = MW-xWatL;
zWatR = zWatL;

% Contour around center (mid) canal
xCanM = MW/2+[   -3.5   -3.0   -2.2   -0.9   -0.0    0.9    2.2    3.0    3.5 ];
zCanM =      [   96.0   88.0   80.0   75.3    74.5   75.3   80.0   88.0   96.0];

xWatM = MW/2+[   -2.2   -0.9   -0.0    0.9    2.2  ];
zWatM =      [   80.0   75.3    74.5   75.3   80.0 ];

milkPoint=   [MW-1.5 D/2 0.5];  % coordates of point where milk was supplied

%% grid, all coordinates in cm

dx  = 1.0; dy=1.8; dz=1.0;
xGr = 0:dx:MW;
yGr = [-D D]/2;
zGr = MH:-dz:0;

gr = gridObj(xGr,yGr,zGr);

%% Zones and plugging of zones into grid

idxMilk = xyzindex(milkPoint,gr);

P = lineGrid([MW/2 D/2 zScrTop; MW/2 D/2 zScrBot],gr);

xWell = [P([1 end end 1]).xm];
zWell = [P([1 end end 1]).zm];

IBOUND = gr.const(0);

IBOUND (inpolyz(gr.XM, gr.ZM, xSand,zSand)  ) = iUnsat;
IBOUND( inpolyz(gr.XM, gr.ZM, xSat  ,zSat ) ) = iSat;
IBOUND( inpolyz(gr.XM, gr.ZM, xCanL ,zCanL) ) = iCanL;
IBOUND( inpolyz(gr.XM, gr.ZM, xCanM ,zCanM) ) = iCanM;
IBOUND( inpolyz(gr.XM, gr.ZM, xCanR ,zCanR) ) = iCanR;
IBOUND( inpolyz(gr.XM, gr.ZM, xWatL ,zWatL) ) = iWatL;
IBOUND( inpolyz(gr.XM, gr.ZM, xWatR ,zWatR) ) = iWatR;
IBOUND( idxMilk )                             = iMilk;
IBOUND( [P.idx] )                             = iWell;


%% Model arrays initially 2D to allow using inpolygon
HK     = gr.const(k);    HK(IBOUND == iUnsat) = k/10; 
VK     = HK;
PEFF   = gr.const(peff);
STRTHD = gr.const(hCanal0);
ICBUND = IBOUND;
STCONC = gr.const(Cfresh);
STCONC(gr.ZM<zIFace)= cSalt;   % Make cells below interface salt

%% Setting boundaries for the stress periods
NPER = getPeriods(basename);

% GHB because head is left free during simulation
zoneValsGHB = { iWatL  'hL' GHBCOND};   % Left  canal GHB

zoneValsCHD = { iWatR  'hR' 'hR';   % Right canal
                iWell  'hW' 'hW';   % Well
                iMilk  'hM' 'hM' };

[GHB, PNTSRC              ] = bcnZone(basename,'GHB',IBOUND,zoneValsGHB,{ Cfresh });

[CHD, PNTSRC(end+(1:NPER))] = bcnZone(basename,'CHD',IBOUND,zoneValsCHD,{ Cfresh; Cfresh; cSalt});

save underneath iCanL iCanR iWatL iWatR iWell iMilk xWell zWell xSand zSand xSat zSat xUnsat zUnsat xCanL zCanL xWatL zWatL xCanR zCanR xWatR zWatR xCanM zCanM xWatM zWatM
