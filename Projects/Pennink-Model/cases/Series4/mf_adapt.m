%% Pennink - density flow in Pennnink's (1915) sand box model
% -- FOURTH SERIES OF EXPERIMENTS (MARCH 1905)
%
%see http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en
%
% in this experiment, Pennink studies simultaneous cFresh and saltwater
% flow from a canal the right and one in the left corner of the sand-box model
% to a canal in the center. The following text is cited from his book, p63:
%
% "B. Phenomena of movement of the liquids of different specific gravity,
% viz water and milk in case the discharge canal is placed in the centre
% and the water symmetrically supplied on both sides"
% Height of the apparatus   = 0.96 m
% Breath of the apparatus   = 0.97 m
% Thickness of the sandmass = 0.018 m
% The reservoir which can be raised and lowered is applied on the right
% side, it was filled with a heavier liquid, viz milk (speicif gravity =
% +/- 1.03), which can enter ton the right hand side from the bottom.
% Water is supplied in the upper corners. Quantity of water = 28 liter
% hourly."
%
% But the surface of the sand has a v-shape with the lowest point in the
% center. THis shape can be used to fix a head gradient at the top of the
% model.
% Initially the water level is above the sand surface and, therefore,
% horizontal, with no flow in the model. The initial level of the milk,
% used as a surrogate for saltwaer, is 25 cm above the bottom of the model.
% The milk is added via a tube near the lower right bottom of the model
% using a reservoir kept at a suitable % elevation to create the initial
% density distribution.
% Next, the water flow is initiated and maintained by adding water to the
% left and right edges on top of the sand surface and extraction it at the sand
% surface in the center of the model. This creates a head gradient equal to
% the inclination of the sand surface, 1:16 according to Pennink on P64.
% Therefore, the flow is driven by the head, which makes the used water
% flow to useless information. It was 28 L/h as mentioned by Pennink.
%
% Because there is no flexibility in this model, there is only one photo in
% Penninks book, showing the situation after some time, which can be
% regarded as the final equilibrium.

% TO 090312 100523 100719

clear variables; close all;

basename='PenninkSeries4';
save name basename

GREP = 'STRESS PERIOD';

fprintf('Tests Pennink (1915), Test series 4, 17 March 1905\n');

%% The grid, the box is 96x96 cm and 1.8 cm thick
cFresh = 0;    % [-] Relative minimum concentration
cSalt  = 1;    % [-] Relative maximum concentration
k      = 86500/(24*60); % [cm/min] calibrated
peff   = 0.38; % [-] effective porosity calibrated
MW     = 97;   % [cm] Width of model. Pennink p6
MH     = 96;   % [cm] Height of model
D      = 1.8;  % [cm] thickness of sand in model

%% Contour of sand body, all coords relative to LL of model sand
xSand =[
   -3.1383
   97.0917
   97.3429
   48.1071
   -0.1238
   -3.1383
];

zSand =[
   -3.4174
   -4.8567
   93.8822
   90.4278
   94.7459
   95.6095
   ];

grad=1./16;  % gradient of sand top acc to Pennink

% top of sand surface, gradient 1/16 acc to Pennink
xTop=   [ 0   0.5 1   ]*MW;
zTop=91+[ 0.5 0   0.5 ]*MW/16;  % 91 cm is top of sand in center of model
    
zIface = 26; % initial elevation of interface

milkPoint = [91.5 1 1.5];  % coordates of point where milk was supplied

%% Grid is always 3D
dx=2.0;              % [m] grid cell width
dy=  D;              % [m] grid cell length = depth of model
dz=2.0;              % [m] grid cell size vertical

xGr=[0:dx:MW MW];    % [m] grid
yGr=[0 dy];          % [m] grid
zGr=[MH:-dz:0 0];    % [m] grid

gr = gridObj(xGr,yGr,zGr);

idxMilk = xyzindex(milkPoint,gr);

iSand = 1;
iMilk = 6;

%% Model arrays

IBOUND = gr.const(0);

IBOUND(inpolyz(gr.XM,gr.ZM,xSand,zSand))  = iSand;

% fixed head at op
IBOUND( diff(IBOUND,1,3) ==1) =-1;

% implant zoneNr into IBOUND
IBOUND( idxMilk) = iMilk;  %milk injection point

HK = gr.const(k);
VK = gr.const(k);

PEFF = gr.const(peff);

STRTHD = XS(ones(gr.Nz,1)*interp1(xTop,zTop,gr.xm));

ICBUND=ones(size(IBOUND));

STCONC = gr.const(cFresh);
STCONC(gr.ZM<zIface) = cSalt;

%% Get period data to set CHD boundaries for MODFLOW and PNTSRC for MTRDMS/SSM
% This is only desired here for the milk injction point which has a fixed
% saltwater head to supply or extrat milk automatically from and to the
% milk reservoir

zoneVals = {iMilk 'hS' 'hS'};

[CHD,PNTSRC] = bcnZone(basename,'CHD',IBOUND,zoneVals,{cSalt});


save underneath xTop zTop
