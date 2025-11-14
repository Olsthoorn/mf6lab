
%% GGOR tool
%

%% What to do

% Show only desired parcels or some overall result
% Same for water budget
% Make a test set and put everything in what needs to be tested
% This should be a separate script that is read in to overwrite the defaults
% Analyze what should be done with tile-drains (dense drain network)
% Analyze what should be done with sub-ditches
% Also check surface runoff and ponding
% See how this is done in MODFLOW, analyze MODFLOW with water table varying
% Document also the testing
% Make paper voor Stromingen met Jos ??
% Python

%%


% Version 2016-01-20
%
% The GGOR tool models the Dutch geohdyrological top system that consists
% of a shallow Holocene aquifer underlain by a
% semi-pervious layer lying on top of a larger regional aquifer.
% The shallow aquiifer is intersected by many dithces in which the surface
% water level is maintained to manage the groundwater system.
%
% We simulate two parallel ditches in this shallow top system underlain by
% the confining bed and the regional aquifer. The ditches may cut through
% the confining layer.
%
% Input to the model is recharge, seepage from below through the confining
% bed and possibly infiltration from the ditches. Extraction consists of
% evapotranspiration and downward leakage plus discharge to the
% ditches.
%
% The specific yield is constant in this model. Making it dependent is postponed
% to a future extension, if required.
%
% When this water table intersects ground surface, surface runoff occurs.
%
% There is no limit to the infiltration capacity of the soil.
%
% Evapotranspiration may be made depdendent on the extincktion depth.
%
% The exchange between ditch and aquifer passes an entry or exit resistance
% that may be different.
%
% Because we enforce that the managed ditch level is the on both sides of
% the modeled parcel, it suffices to simuulate only half parcel.
%
%
% TO 100905 101021 151112 160120

%clear variables; close all;

basename='Noo_GGOR_2015_opgeknipt';  % 

%% Read meteo time series

TEST = false;

tne    = getMeteo('../Meteo/PE-00-08.txt');

if TEST == true
    tne(:,3) = 0;

    tt =  365*[0  1  2  3  4  5  6 7];
    pp = 0.01 *[0  1 -1  1 -1  1 -1 -1];

    for it=1:numel(tt)
        tne((180+tt(it)):end,2) = pp(it);
    end

end

%% set wat winter is
DV     = datevec(tne(:,1));        % tne(:,1) is simulation time
winter = DV(:,2)>=10 | DV(:,2)<=3; % winter is a logical vector

%% Maximum parcel width. Parcles > Lmax have subdrainage systems
Bmax = 50; % default 50 m,  [m] parcels > Bmax must have intermediate ditches that are unmapped but visible on areal photos

%% Read the data for all parcels
% dbfRead yields the columns in the dbfFile and the
% fieldnames corresponding to these columns

dbfCol  = dbfread(['../GGORdata/' basename]);

fldNms  = deblank({dbfCol.fieldname}); %remove blanks from field names

%% Cleanup data to work with

% Rule to Fetch the dbfCol
Fetch  = @(fld) dbfCol(strmatchi(fld,fldNms,'exact')).values;

% Rule to convert a scalar to a vector for all parcels
parcels  = 1:numel(dbfCol(1).values); % set of parcels to consider
Nparcels = numel(parcels);
const  = @(val) val * ones(Nparcels,1);

%% Remove all parcels from the database with incomplete dbfCol (having NaNs)

Iparcels = true(numel(dbfCol(1).values),1);

% make sure all columns have values
for j=1:numel(fldNms)
    Iparcels = Iparcels & ~isempty(dbfCol(j).values);
end

hk1 = Fetch('K');

D1     = Fetch('DIKTE_DEK');   % thickness of holocene top aquifer (cover layer)
c      = Fetch('C_DEK');   % [ d ] parcel resistance
C_DEK_MIN=1;       % default minimum resistance cover layer 
D1       = max(0.1,D1);
C        = max(C_DEK_MIN,c);
vk1      = 0.5 * D1./C; 
sy1       = Fetch('MU');         % [ - ] Storage coefficient layer 1

% %% Get kh1 and sy1 according to their BOFEK codes
% BOFEK    = Fetch('BOFEK');
% [bofekHdr,bofek] = getExcelData(basename,'BOFEK','HOR');
% 
% ibo = strmatchi('BOFEK',bofekHdr);
% ikh = strmatchi('kh',bofekHdr);
% isy = strmatchi('sy',bofekHdr);
% 
% % one should check that all parcels have obtained a value
% hk1 = const(NaN);
% vk1 = const(NaN);
% sy1 = const(NaN);
% for j=1:size(bofek,1)
%     I = BOFEK == bofek(j,ibo);
%     if ~isempty(I)
%         hk1(I) = bofek(j,ikh);
%         vk1(I) = bofek(j,ikh); %JB: verticale anisotropie toevoegen?
%         sy1(I) = bofek(j,isy);
%     end
% end
% if any(isnan(hk1))
%     error(['Not all parcels have been attributed a hk1 and sy1 value\n',...
%        'check that the BOFEK sheet in excel (basename) is complete, in that\n',...
%        'all used BOFEK numbers have their corresponding hk and sy values specified.']);
% end

%% Percel width (formule JW Voort, notitie 08.002684)
AREA = Fetch('AREA');      % [m2] true parcel area from GIS
O = Fetch('perimeter'); % [m2] true parcel perimeter from GIS
D = O.^2-16.*AREA;         % discriminant
I = D>=0;               % discriminant>0 ? --> real solution
B = NaN(size(I));       % parcel width
L = NaN(size(I));       % parcel length
B( I)  = (O(I) - sqrt(D(I)))/4;  % width, smallest of the two values
L( I)  = (O(I) + sqrt(D(I)))/4;  % length, largest of the two values
B(~I)  = sqrt(AREA(~I));   % if no real solution --> assume square
L(~I)  = sqrt(AREA(~I));   % same, for both width and length

B(B>Bmax) = Bmax;
b      = B/2;           % half parcel width
bmin     = 10; % min width of parcel
Iparcels = Iparcels & b>bmin;
Iparcels = Iparcels & hk1 <1000;
ZP     = Fetch('ZP_2015'); % summer ditch water level 
Iparcels = Iparcels & ZP<0;

landuse = Fetch('TYPELANDGE');%landgebruikskaart Waternet
usezones=['overig','bebouwd gebied'];
Iparcels=Iparcels & ~ismember(landuse,usezones);

% FID    = Fetch('FID4');% beperk het aantal percelen
% Iparcels=Iparcels & FID<50;% beperk het aantal percelen


%% Selecteer parcels with larger than mininum width
% ZONE = Fetch('ZONE_CODE');

% useZones = [3, 4, 5, 6, 11, 12, 13, 23, 32, 38, 44, 52]; % in polders BKP en UHP

% Iparcels = Iparcels & ismember(ZONE,useZones); 
% landuse = Fetch('WN_LGN6');%landgebruikskaart Waternet
% 
% useZones = [2,3,5]; %akkerbouw, tuinbouw, gras
% Iparcels = Iparcels & ismember(landuse,useZones); 
% 
Iparcels = Iparcels & hk1>0;
b = b(Iparcels);
O = O(Iparcels);
AREA = AREA(Iparcels);
hk1=hk1(Iparcels);
vk1=vk1(Iparcels);
sy1=sy1(Iparcels);
B=B(Iparcels);
C=C(Iparcels);
D=D(Iparcels);
D1=D1(Iparcels);
L=L(Iparcels);
ZP=ZP(Iparcels);
c=c(Iparcels);
%% If you want to further select a subset of the parcels then do that here
Nsubset = 10;        
if TEST==true
    I = find(Iparcels,Nsubset,'first');
    Iparcels(I(end)+1:end) = false;
end    

%% Only keep the fields with non NaN numeric data
for ic = 1:numel(fldNms)
    dbfCol(ic).values = dbfCol(ic).values(Iparcels);
end

%% Redefine rule to Fetch the dbfCol, because dbfCol has changed
Fetch  = @(fld) dbfCol(strmatchi(fld,fldNms,'exact')).values;
Nparcel = sum(Iparcels);


%% Rule to convert a scalar to a vector for all parcels

const  = @(val) val * ones(Nparcel,1);

%% Choose a subset of the parcels base for exercising with the GGOR tool

FID    = Fetch('FID4');
%AREA   = Fetch('AREA');    % [m2] true parcel area from GIS
%O      = Fetch('perimeter'); % [m2] true parcel perimeter from GIS

xCtr   = Fetch('X');     % parcel center x
yCtr   = Fetch('Y');     % parcel center y
%ZP     = Fetch('ZP_2015'); % summer ditch water level 
WP     = Fetch('WP_2015'); % winter ditch water level
GP     =(ZP+WP)/2;
%D1     = Fetch('DDEK');   % thickness of holocene top aquifer (cover layer)
DCB    = const(0.01);     % dummy confining bed (CB) thickness below cover layer (represents its vertical resistance)
D2     = const(30);       % thickness of pleistocene aquifer (2dn layer)
%BOFEK  = Fetch('BOFEK');  %code voor bodemtype BOFEK bodemkaart.

%hk1    = Fetch('hk1');     % [m/d] kh in holocene layer
%sy1    = Fetch('MU2');   % parcel Sy

hk2    = const(30);       %[m/d] kh default in regional aquifer
%c      = Fetch('Cdek');   % [ d ] parcel resistance
vk2    = hk2;         % theo: hk2/500 is dit een test?? [ m/d ] default vertical cond regional aquifer
                        % bedenk dat de eerste meters vrijwel altijd Boxtel
                        % zand is met een k-waarde van 1-2 m/dag en een
                        % verticale anistropie van 5 of hoger...
sy2    = const(0.22);     % in case regional aquifer becomes phreatic
ss1    = const(1e-5);     % [1/m] elastic storativity --> optie gebruik dit voor waterberging op maaiveld
ss2    =1e-5* D2;     % [1/m] same regional aquifer
q      = Fetch('Qflexkwar1'); % [m/d] upward seepage
phi    = Fetch('PHI2');       % [ m ] head in regional aquifer (not used)
AHN    = Fetch('AHN2010MED'); AHN(abs(AHN)>10)=0; % ground elevations

dw     = const(2.0);     % [ m ] default width of ditch
drnFac = const(2.0);     %  cIn/cEx;
cEx    = const(2.0);     % [ d ] exfiltration resistance ditch
cIn    = cEx .* drnFac;  % [ d ] default entry resistance of ditch
dd     = const(0.6);     % [ m ] default ditch depth = 0.6 m
omega  = dw+2*dd;        % [ m ] default circumference of ditch--> check of dit goed gaat verderop
cdr    = const(0.1);     % [ d ] default drain resistance (as surface resistance)
bodcod1= Fetch('BODCOD1'); % de bodemtypering volgens Alterra.

%% Test
if TEST == true
    b = round(b);
    L = 2*b;
    %dd = 0.5*(WP+ZP) - (AHN-(D1+1)); %ditches through cover-layer
    q = q*10;
end
%% compute Omeg1 and Omeg2
delta    = (ZP+WP)/2 - dd - (AHN - D1); % remaining thickness cover layer below ditch bottom

% ditch width
dw1      = 2 * min(0.5*dw,max(0,delta));
dw2      = dw - dw1;

% ditch depth
%dd2      = max(0,-delta); %Theo O
dd2      =min(dd,max(0,-delta)); %Josbe een maximum ingebouwd omdat delta geen limiet heeft met de diepte.
dd1=max(0.1,dd-dd2);% josbe een minimum ingebouw om te voorkomen dat geen drainage optreedt in laag 1.
%dd1      = dd - dd2;% Theo O

% wet circumference
Omeg1    = 2 * dd1 + dw1;
Omeg2    = 2 * dd2 + dw2;

%% Conductivties and lambda?
kD1      = hk1 .* D1;
lambda   = sqrt(kD1.*c);
k1z      = D1./c;
k1       = sqrt(hk1.*vk1);
k2       = sqrt(hk2.*vk2);


%% Specify ditch levels over time
hDitch = WP * winter' + ZP * ~winter';   % set ditch level

%%  Coordinates for the model grid
dx  = 2;  % cell width choice, default 1 m

% cell boundary coordinates using a tiny left cell to set boundary
% conditions at virtually x = 0. (May be altered layer).
xGr = [ -0.001 0:dx:round(max(b)) ];

% One parcel in each row of the model (set anisotropy to 10-10 in LAY sheet) 
yGr = 0:1:Nparcel;

Nx  = numel(xGr)-1;
Ny  = numel(yGr)-1;

% Rule to convert a parcel vector every model cell in one layer
lay = @(var) var * ones(1,Nx);

%% Z of model grid

deltaZ = 0.001;  % use this for the confining bed at bottom of 1st aquifer

Z   = NaN(Ny,Nx, 3); % allocate memory for Z grid

Z(:,:,1) = lay(AHN);                     % ground surface elevation
Z(:,:,2) = Z(:,:,1) - lay(D1);           % bottom of cover layer
Z(:,:,3) = Z(:,:,1) - lay(D1) - deltaZ;  % bottom of confining bet (=cover layer)
Z(:,:,4) = Z(:,:,1) - lay(D1) - lay(D2); % bottom of regional aquifer

%% Generate grid Object

% guarantee MINDZ not smaller than thickness of (dummy) confining bed
% use LAYCBD to define confining bed at bottom of cover layer
gr = gridObj(xGr,yGr,Z,'MINDZ',deltaZ,'LAYCBD',[1 0]);

%% Id's for set and remember cell properties
iDRN  = 1;    % zone nr of drains in the model
iGHB1 = 2;    % zone nr of ditchs in the model
iGHB2 = 4;    % same in regional aquifer
iRIV1  = iGHB1; % zone nr of river cells in the model (same as iGHB)
iRIV2  = iGHB2; % zone nr of river cells in the model (same as iGHB)
iWVP2 = 3;    % zone nr of regional aquifer (where to inject the seepage)
iWEL  = iWVP2;% zone nr where to inject the seepage (same as regional aquifer)

%% Build IBOUND array using zone numbers to easily identify zones
IBOUND = gr.const(1);
IBOUND(:,    :,1) = iDRN;    % zone with drains (ground surface)
IBOUND(:,    :,2) = iWVP2;   % second aquifer
IBOUND(:,    1,1) = iGHB1;   % location of ditch, overwrite iDRN
IBOUND(:,    1,2) = iGHB2;   % location of ditch, overwrite iWVP2

% iWEL and iRIV not necessary to set in IBOUND because equal to iWVP2 and iGHB

%% Set active width of each model equal to parcel width b
% This is done by making cells beyond b inactive (IBOUND==0)

for iy=1:gr.Ny
    IBOUND(iy,gr.xm>b(iy),:)=0;
end

% mark which cells are active
active = logical(IBOUND(:,:,1)>0);

%% Build model arrays
HK            = gr.const(NaN);
HK(:,:,1)     = lay(hk1);
HK(:,:,2)     = lay(hk2);

VK            = gr.const(NaN);
VK(:,:,1)     = lay(vk1);
VK(:,:,2)     = lay(vk2);

VKCB          = lay(deltaZ./c);

SY            = gr.const(NaN);
SY(:,:,1)     = lay(sy1);
SY(:,:,2)     = lay(sy2);

SS            = gr.const(NaN);
SS(:,:,1)     = lay(ss1);
SS(:,:,2)     = lay(ss2);

% Notice if laycon == 0 (e.g. to verify GGOR)
% make sure Sy is used, put it in Ss of first layer
[layHdr,layVals] = getExcelData(basename,'LAY','horizontal');
if layVals(1,strmatchi('LAYCON',layHdr))==0
    SS(:,:,1)    = SY(:,:,1)./gr.DZ(:,:,1);
end

%% Initial heads
STRTHD        = gr.const(NaN);
STRTHD(:,:,1) = lay(hDitch(:,1));
STRTHD(:,:,2) = lay(phi);

%% Stress periods to generate RECH and EVTR for MODFLOW
% Use mf2005 because mf2k has a limit of 1000 stress periods

[~,~,NPER]=getPeriods(basename);

% Check that NPER matches length of time series in tne
if size(tne,1) ~= NPER
    error('tne length = %d does not mathc NPER = %d\n',size(tne,1),NPER);
end

%% Recharge MODFLOW specfied here (make sure INRCH == 0 in PER)

% The recharge input is a Ny*Nx*Nper array but can be one of size (1,1,NPER)
% in which case each value is the recharge for the whole model
% for the corresponding stress period.
RECH   = ones(1,1,NPER);
RECH(:)= tne(:,2);

%% Evapotranspiration for MODFLOW specified here
%  EVTR separate from RECH allows reduction of evapotr. in dry spells
EVTR   = ones(1,1,NPER);
EVTR(:)= tne(:,3);

% SURF   = repmat( {AHN - 30.0} , NPER ,1);  % Zie PER sheet
% Make sure the firs inSurf == -2 so that SURF in the PER sheet
% is interpreted as the distance below the top of the model.

%% Drains at ground surface to compute surface runoff
Cdrn   = gr.AREA ./ lay(cdr);   % Drain conductance
Idrn   = find(IBOUND == iDRN);  % Drain locatioins in network
LRCdrn = cellIndices(Idrn,gr.size,'LRC');

%% General Head Boundaries (GHB) in the left-most cell to simulate the ditch
% We need GHB to include entry and exit resistances separately with the ditches; 
% we ignore radial resistance for the time being %Josbe: deze opmerking kan weg, radial res. staat hieronder:
rho1    = max(0,2./pi ./ k1 .* log(D1 ./ Omeg1)); 
rho2    = max(0,2./pi ./ k2 .* log(D2 ./ Omeg2)); 

Rex1    = cEx ./ (0.5*Omeg1);
Rex2    = (cEx + max(0,delta./k1z)) ./ (0.5*Omeg2);

Rghb1   = cIn ./ (0.5*Omeg1);
Rghb2   = (cIn + max(0,delta./k1z)) ./ (0.5*Omeg2);

Rriv1   = (Rex1+rho1).*(Rghb1+rho1)./(Rghb1-Rex1)-rho1; % Josbe: p23
Rriv2   = (Rex2+rho2).*(Rghb2+rho2)./(Rghb2-Rex2)-rho2; % Josbe: p23

Cghb1   = gr.dy ./ (Rghb1 + rho1); % GHB conductance <-- entry resistance
Cghb2   = gr.dy ./ (Rghb2 + rho2); % GHB conductance <-- entry resistance
Cghb2(Omeg2==0) = 0;

Ighb1   = find(IBOUND==iGHB1);  % global indices of ditch
Ighb2   = find(IBOUND==iGHB2);  % global indices of ditch

LRCghb1 = cellIndices(Ighb1,gr.size,'LRC');
LRCghb2 = cellIndices(Ighb2,gr.size,'LRC');

%% RIV to simulate the exit resistance in combination with GHB
Criv1   = gr.dy ./ ( Rriv1 + rho1);
Criv2   = gr.dy ./ ( Rriv2 + rho2);
Criv2(Omeg2==0) = 0;

Iriv1   = find(IBOUND==iRIV1);
Iriv2   = find(IBOUND==iRIV2);

LRCriv1 = cellIndices(Iriv1,gr.size,'LRC');
LRCriv2 = cellIndices(Iriv2,gr.size,'LRC');

%% Wells to simulate seepage from the bottom aquifer or vice versa
qWEL        = gr.const(0);  % need two layer because well is in layer 2
qWEL(:,:,2) = gr.AREA .* lay(q); % set prescribed seepage in layer 2
Iwel        = find(IBOUND==iWEL);     % Global indices of well
LRCwel      = cellIndices(Iwel,gr.size,'LRC');

%% need a column of ones of given size to multiply scalars with
uDRN  = ones(size(Idrn));
uWEL  = ones(size(Iwel));
uGHB1 = ones(size(Ighb1));
uGHB2 = ones(size(Ighb2));
uRIV1 = ones(size(Iriv1));
uRIV2 = ones(size(Iriv2));

%% Specify boundary conditions (stresses) for the MODFLOW model

for iPer = NPER:-1:1
    if iPer == 1
        DRN{iPer,1} = [iPer*uDRN LRCdrn Z(Idrn) Cdrn(Idrn)];
        WEL{iPer,1} = [iPer*uWEL LRCwel qWEL(Iwel)];
    else
        DRN{iPer,1} = [-iPer ones(1,5)]; % as previous
        WEL{iPer,1} = [-iPer ones(1,4)]; % as previous
    end
        
    if iPer==1 || (winter(iPer) * winter(iPer-1) == false) % change of season
        if winter
            GHB{iPer,1} = [[iPer*uGHB1 LRCghb1 WP Cghb1]   ;...  % first aquifer
                         [iPer*uGHB2 LRCghb2 WP Cghb2]];       % second aquifer
            RIV{iPer,1} = [[iPer*uRIV1 LRCriv1 WP Criv1 WP];...  % first aquifer
                         [iPer*uRIV2 LRCriv2 WP Criv2 WP]];    % second aquifer
        else
            GHB{iPer,1} = [[iPer*uGHB1 LRCghb1 ZP Cghb1]   ;[iPer*uGHB2 LRCghb2 ZP Cghb2]];
            RIV{iPer,1} = [[iPer*uRIV1 LRCriv1 ZP Criv1 ZP];[iPer*uRIV2 LRCriv2 ZP Criv2 ZP]];
        end
        % remove lines with zero conductance, they only take file space but
        % are no use for the model
        GHB{iPer,1} = GHB{iPer,1}(GHB{iPer,1}(:,end  )~=0,:);
        RIV{iPer,1} = RIV{iPer,1}(RIV{iPer,1}(:,end-1)~=0,:);
    else
        GHB{iPer,1} = [-iPer ones(1,5)]; % as previous
        RIV{iPer,1} = [-iPer ones(1,6)]; % as previous
    end    
end

%% Simulation using analytical solution with constant layer thickness and prescibed flux

% Be aware though, that this simulation totally ignores surface runoff
% and evaporation reduction durign dry spells.
% Therefore the numerical and analytical must differ unless the numerical
% settings are set such that the analytical circumstances are guaranteed
% in the numerical model.
% Lastly, the analytical solution works with constant aquifer thickness (to
% avoid zero thickness in case of large downward flow.


% resistance cover layer %josbe: and radial resistance; log() met ln() zijn.
wEx     = 2 * cEx .* D1 ./ Omeg1 + max(0,2/pi .* D1./k1 .* log(D1 ./ Omeg1)); % exit resistance (as a wall of height D1) 
wIn     = 2 * cIn .* D1 ./ Omeg1 + max(0,2/pi .* D1./k1 .* log(D1 ./ Omeg1)); % entry resistance(as a wall of height D1)

% reistance ditches in regional aquifer

gammaEx     = (cEx + max(0,delta)./k1z) .* B./Omeg2 + B./(pi*k2) .* log(D2./Omeg2);
gammaIn     = (cIn + max(0,delta)./k1z) .* B./Omeg2 + B./(pi*k2) .* log(D2./Omeg2);

cGammaEx = c .* ( b./c .* wEx./D1 + b./lambda .* coth(b./lambda) - 1);
cGammaIn = c .* ( b./c .* wIn./D1 + b./lambda .* coth(b./lambda) - 1);

cGcgEx   = cGammaEx + c + gammaEx;
cGcgIn   = cGammaIn + c + gammaIn;

Tex      = sy1 .* cGammaEx;  % Time constant when exfiltrating
Tin      = sy1 .* cGammaIn;  % Time constant when infiltrating

%% Analytical simulation

% Allocate space to store the analyical heads (average head in cross sections)
hAnalytic      = NaN(Nparcel,NPER);  % head in the shallow aquifer
hAnalytic(:,1) = hDitch(:,1);        % initialize first values
q1             = NaN(Nparcel,NPER);  % seepage through cover layer
q2             = NaN(Nparcel,NPER);  % seepage through cutting ditches
phi            = NaN(Nparcel,NPER);  % seepage through cutting ditches
eEx            = NaN(Nparcel,1);
eIn            = NaN(Nparcel,1);
N              = tne(:,2)-tne(:,3);  % net recharge
% Simulate all parcels simultaneously while differentiating between in- and exfiltration

Dt       = diff(tne(:,1));

for it=1:length(Dt);
    
    Iex = hAnalytic(:,it)>hDitch(:,it); % logial indicating which parcels exfiltrate at time step it
    
    q1( Iex,it) = (q( Iex) .* gammaEx( Iex) - N(it) .* cGammaEx( Iex)) ./ cGcgEx( Iex);
    q1(~Iex,it) = (q(~Iex) .* gammaIn(~Iex) - N(it) .* cGammaIn(~Iex)) ./ cGcgIn(~Iex);
    
    q1( Omeg2 == 0, it ) = q( Omeg2 == 0);
    
    q2(:,it) =  q(:) - q1(:,it);
    
    eEx(:) = exp(-Dt(it)./Tex);
    eIn(:) = exp(-Dt(it)./Tin);
    
    hAnalytic(Iex,it+1) = hDitch(Iex,it) + ...
            ( hAnalytic( Iex,it) - hDitch( Iex,it) ) .* eEx( Iex) +  ...
             cGammaEx( Iex) .* ( N(it) + q1( Iex,it) ) .* ( 1 - eEx( Iex) ); 
         
    hAnalytic(~Iex,it+1) = hDitch(~Iex,it) + ...
            ( hAnalytic(~Iex,it) - hDitch(~Iex,it) ) .* eIn(~Iex) + ...
             cGammaIn(~Iex) .* ( N(it) + q1(~Iex,it) ) .* ( 1 - eIn(~Iex) );
         
    % verification equation phi-hLR in Lyx text
    phi(Iex,it)  = hDitch( Iex,it) + N(it) * cGammaEx( Iex) + q1( Iex,it).*(cGammaEx( Iex) + c( Iex));
    phi(~Iex,it) = hDitch(~Iex,it) + N(it) * cGammaIn(~Iex) + q1(~Iex,it).*(cGammaIn(~Iex) + c(~Iex)) ;


    fprintf('.'); if rem(it,100)==0, fprintf('\n'); end % show your're busy
end
fprintf('\n');

%% Analytically compute the final head for steady state 
% This is for checking GGOR, may be skipped here, see steadySolution.m

% Ditch resistance in last time step
w = NaN(Nparcel,1);
w( Iex) = wEx( Iex);
w(~Iex) = wIn(~Iex);

b_lam = b./lambda;

% x-coordinate along the analytical cross section
x  = bsxfun(@minus, b, gr.xm);

% Analytical steady state solution [m]
hx = NaN(Nparcel,gr.Nx);

hx( Iex, :) = bsxfun(@times, phi( Iex,end) + N(end).*c( Iex), active( Iex,:)) - ...
              bsxfun(@times, (phi( Iex,end) + N(end).*c( Iex) - hDitch( Iex,end ) ) ./ ...
         ((lambda( Iex)./c( Iex) .* w( Iex)./D1( Iex)).*sinh(b_lam( Iex)) + cosh(b_lam( Iex))), ...
          cosh( bsxfun(@times, x( Iex,:), 1./lambda( Iex)) ) );
hx(~Iex, :) = bsxfun(@times,  phi(~Iex,end) + N(end).*c(~Iex), active(~Iex,:)) - ...
              bsxfun(@times, (phi(~Iex,end) + N(end).*c(~Iex) - hDitch(~Iex,end)) ./ ...
         ((lambda(~Iex)./c(~Iex) .* w(~Iex)./D1(~Iex)).*sinh(b_lam(~Iex)) + cosh(b_lam(~Iex))), ...
          cosh( bsxfun(@times, x(~Iex,:), 1./lambda(~Iex)) ) );

% Analytical steady state solution for flow to or from the ditch [m2/d]
Q = NaN(Nparcel,1);
Q( Iex) = (phi( Iex,end) + N(end).*c( Iex) - hDitch( Iex,end)) .* sinh(b_lam( Iex)) ./ ...
    ((w( Iex)./D1( Iex)).*sinh(b_lam( Iex)) + (lambda( Iex)./kD1( Iex)) .* cosh(b_lam( Iex)));
Q(~Iex) = (phi(~Iex,end) + N(end).*c(~Iex) - hDitch(~Iex,end)) .* sinh(b_lam(~Iex)) ./ ...
    ((w(~Iex)./D1(~Iex)).*sinh(b_lam(~Iex)) + (lambda(~Iex)./kD1(~Iex)) .* cosh(b_lam(~Iex)));

%% Save whatever non-modflow data is necessary in mf_analyze
save underneath FID sy1 AHN tne phi q kD1 hk1 vk1 D1 w c b omega lambda hDitch hAnalytic AREA L bodcod1 GP;
