%% Analyzing the simulation result of the shallow Dtuch top system.
%
% TO 100823 151116 151126

clear variables; close all;

load name
load(basename)
load underneath

m2mm = 1000; % convert from m to mm
ha   =10000; % hectare in m

fsz      = 12;
defaults = {'nextPlot','add','xGrid','on','yGrid','on','fontsize',fsz};

active = logical(IBOUND(:,:,1));

%% Get heads of first layer and second layer separately

top = 1;
bot = 2;

H    =readDat([basename '.HDS'],'','',top); % Get only the first layer of heads
Phi = readDat([basename '.HDS'],'','',bot);  % Get head of reagional aquifer

Nt   = size(tne,1);
Nsec = gr.Ny;

% set H at inactive cells to zero (prevent NaNs)
for it=1:Nt, H(it).values(~active) = 0; end

%% Compte the average numeric head over cross sections

Area     = sum(gr.AREA .* active,2); % only what's active !

% Make room and compute cross section average heads (which is what we want for GGOR)
hNumeric = zeros(Nsec,Nt);
for it=1:Nt
    hNumeric(:,it) = sum(H(it).values .* gr.AREA, 2)./Area;
end

%% Compute GxG separately for numerical and analytical model

%[GLGn,GVGn,GHGn] = getGXG(hNumeric ,tne(:,1),'plot',[1 2],'title','numeric'); % numerical
%[GLGa,GVGa,GHGa] = getGXG(hAnalytic,tne(:,1),'plot',[1 2],'title','analytic'); % analytic
[GLGn,GVGn,GHGn] = getGXG(hNumeric ,tne(:,1),'plot',[1 2],'title','numeric'); % numerical
[GLGa,GVGa,GHGa] = getGXG(hAnalytic,tne(:,1),'plot',[1 2],'title','analytic'); % analytic
%% Compute area-weighted GXG over all parcels in the model

% Notice that AREA (all capital) is the vector with parcel areas from the database

GLGn_avg = sum(GLGn .* AREA,1) / sum(AREA);
GVGn_avg = sum(GVGn .* AREA,1) / sum(AREA);
GHGn_avg = sum(GHGn .* AREA,1) / sum(AREA);

GLGa_avg = sum(GLGa .* AREA,1) / sum(AREA);
GVGa_avg = sum(GVGa .* AREA,1) / sum(AREA);
GHGa_avg = sum(GHGa .* AREA,1) / sum(AREA);

%% plot GxG for all parcels from the analytic and the numeric model

fsz = 12;
figure('position',screenPos(0.8));
ax = axes(defaults{:});
xlabel(ax,'Cross section number');
ylabel(ax,'head (above datum)');
title(ax,'GxG of the computed cross sections (GGOR tool, MODFLOW and analytical)');

plot(ax,1:Nsec,GLGn,'r-');
plot(ax,1:Nsec,GVGn,'g-');
plot(ax,1:Nsec,GHGn,'b-');
plot(ax,1:Nsec,GLGa,'r:');
plot(ax,1:Nsec,GVGa,'g:');
plot(ax,1:Nsec,GHGa,'b:');
legend(ax,'GLGn','GVGn','GHGn',...
          'GLGa','GVGa','GHGa');
    
%% Compare analytic and numeic cross-section averaged heads versus time
% figure('position',screenPos(0.8));
% 
% ax1 = axes(defaults{:});
% 
% xlabel(ax1,'time');
% ylabel(ax1,'head (above datum)');
% title(ax1,{'Comparison of analytically and numerically computed heads';...
%           'for selected sections; analytic blue, numeric red'});
% 
% Isections = [1 2 3]; % sections to show here
% 
% plot(ax1,tne(:,1),hAnalytic(Isections,:),'b');
% plot(ax1,tne(:,1),hNumeric( Isections,:),'r');

% datetick;

return

%% === water budgets ===== water budgets ===== water budgets =====

% Use these two rules for an intermediate pause in case budget computation
% takes very long. If no budgets are needed the folowing section may be
% skipped by placing a return statement here:
%return

%fprintf('press RETURN to continue\n');
%pause;

LABELS = {'EVT'  'ET'            'EVTR'   'QET>0'  'EVTR'   'QET<0'  'y';
          'RCH'  'RECHARGE'      'RCH'    'QRCH>0' 'RCH'    'QRCH<0' 'm';
          'STO'  'STORAGE'       'STO'    'QSTO>0' 'STO'    'QSTO<0' 'b';
          'CHD'  'CONSTANTHEAD'  'CHD'    'QCHD>0' 'CHD'    'QCHD<0' 'w';
          'WEL'  'WELLS'         'UPWSPG' 'QWEL>0' 'DWNSPG' 'QWEL<0' 'g';
          'DRN'  'DRAINS'        'RUNOFF' 'QDRN>0' 'RUNOFF' 'QDRN<0' 'w';
          'RIV'  'RIVERLEAKAGE'  'RIV'    'QRIV>0' 'RIV'    'QRIV<0' 'r';
          'GHB'  'HEADDEPBOUNDS' 'GHB'    'QGHB>0' 'GHB'    'QGHB<0' 'c';
          'LKG'  'FLOWLOWERFACE' 'LKG'    'QLOW>0' 'LKG'    'QLOW<-' 'k'}; % seepage cover laye3

%% Read cell by cell flows form budget files

% Use LABELS(:,2) to save some space in memory, only read these data
% Notice however, that all necessary data must be read
B = readBud([basename '.BGT'],LABELS(:,2)');

%% Side step:  compare the analytic and numeric steady state solution
% This is mainly for testing the GGOR tool
steadySolution;  % testing is done in separate script to easily skip it

% If you one to continue with the budgets, outcomment the return statement
%return
  
%% Narrative water budgets
% For each stress period the budget file hold a 3D array for each cell by
% cell flow, an amazing amount of data.
% 
% First we sum the CBC data over every cross section, so that
% we end up with the total for the DRN, GHB ect flow for every cross section
% and stress period.
%
% We compute the water budget for every cross section over time from this.
%
% Then we can combine the water budget of arbitrary sets of cross sections.
% We may for instance compute the overall running water balance for
% the entire AREA by summing over all cross sections weighted by
% the actual parcel AREA given in the database. We could do so for every
% CBC flow separately, to find the total runoff (= total DRN), the total
% recharge, the total ditch discharge (= GHB) and the total seepage (= given WEL).
% By adding the index of the desired subareas in the database, this
% totaling can be done automatically.

%% Terms for the water budget see LABELS
% The terms of interest in this model are
% RCH  -- recharge
% EVT  -- evapotranspiration
% GHB  -- flow from the ditches
% RIV  -- exfiltration along with GHB
% The exchange between ground and ditch is just GHB+RIV
% DRN  -- overflow (surface runoff)
% WEL  -- fixed seepage flow injected in second aquifer
% LKG  -- flow lowerface, direct exchange between cover layer and regional aquifer
% STO  == storage during time step
%
% To readily compare these values, we may compute the values
% directly per unit of ground surface [L/T] = [L^3/T/L^2]

%% Allocate memory to store budget terms.

% NSec = cross section, Nt is time
stress = ['STO' whichStresses(LABELS(:,1))' 'LKG'];

for i = 1:numel(stress)
    Q1.(stress{i}) = zeros(Nsec,Nt);
    Q2.(stress{i}) = zeros(Nsec,Nt);
end

% stress position (row) in LABELS)
Istress = zeros(size(stress));
for i=1:numel(stress)
    Istress(i) = strmatchi(stress{i},LABELS(:,1));
end

%% Fill the arrays using the flow terms from the budget file in [L/T] units

% We must realize that not all labels must be present at all time steps.
% However they are in the GGOR tool and they correspond to the stresses in
% stress.

% Sum over columns and layers to get totals for all cross sections as a
% funcion of time

% rule to sum cell by cell flow of given type
add1 = @(it,type) sum(B(it).term{type}(:,:,  1),2) ./ Area;
add2 = @(it,type) sum(B(it).term{type}(:,:,end),2) ./ Area;

for it=1:Nt
    for i=1:numel(Istress)
        iterm = strmatchi(LABELS{Istress(i),2},B(it).label);
        Q1.(stress{i})(:,it) = add1(it,iterm);
        Q2.(stress{i})(:,it) = add2(it,iterm);
        
        % LKG of first layer is - flow FTF of second layer
        if strcmpi(stress{i},'LKG')
            Q2.(stress{i})(:,it) =  Q1.(stress{i})(:,it);
            Q1.(stress{i})(:,it) = -Q1.(stress{i})(:,it);
        end
    end
end

%% Discharges per unit ground surface [L/T] over the entire GGOR region

% The relative contribution of parcel to the stress on the stress value over
% all values would be its stress per unit area in hte model multiplied by
% its real world area of the parcel divided by the sum of the area of all
% parcels
f = AREA ./ sum(AREA) * ones(1, Nt);

for i = 1:numel(stress)
    Qtot1.(stress{i}) = sum(Q1.(stress{i}).* f, 1);
    Qtot2.(stress{i}) = sum(Q2.(stress{i}).* f, 1);
end

%% Total these flows to see if the water budget adds op to zero (or almost)

QTOT1 = zeros(size(Qtot1.(stress{1})));
QABS1 = zeros(size(Qtot1.(stress{1})));
QTOT2 = zeros(size(Qtot2.(stress{1})));
QABS2 = zeros(size(Qtot2.(stress{1})));

for i=1:numel(stress)
    QTOT1(:,:) = QTOT1 +           Qtot1.(stress{i});
    QABS1(:,:) = QABS1 + 0.5 * abs(Qtot1.(stress{i}));
    QTOT2(:,:) = QTOT2 +           Qtot2.(stress{i});
    QABS2(:,:) = QABS2 + 0.5 * abs(Qtot2.(stress{i}));
end
%% Display these flows summed over all times for the cross sections
% This is to show that the water budget matches for all cross section of all times

fprintf('%10s',stress{:},'[mm/d]'); fprintf('\n');

A1 = zeros(Nsec,numel(stress));
A2 = zeros(Nsec,numel(stress));

for i=1:numel(stress)
    A1(:,i) = mean(Q1.(stress{i}),2);
    A2(:,i) = mean(Q2.(stress{i}),2);
end
display(m2mm * A1);
fprintf('\n');
display(m2mm * A2);

%% Show running water budget summed over all AREA-weigted cross sections 
%  selected from the database

QPOS1  = zeros(numel(stress),Nt);
QNEG1  = zeros(numel(stress),Nt);
QPOS2  = zeros(numel(stress),Nt);
QNEG2  = zeros(numel(stress),Nt);

for i=1:numel(stress)
   QPOS1(i,:) = Qtot1.(stress{i});
   QNEG1(i,:) = Qtot1.(stress{i});
   QPOS2(i,:) = Qtot2.(stress{i});
   QNEG2(i,:) = Qtot2.(stress{i});
end

QPOS1(QPOS1<0) = 0;
QNEG1(QNEG1>0) = 0;
QPOS2(QPOS2<0) = 0;
QNEG2(QNEG2>0) = 0;
   
%%
figure('position',screenPos(0.8));

ax1 =subplot(2,1,1,defaults{:},'xlim',tne([1 end],1));
ax2 =subplot(2,1,2,defaults{:},'xlim',tne([1 end],1));

xlabel(ax1,'tijd'); ylabel(ax1,'m/d');
title( ax1,'Voortschrijdende waterbalans delaag');
xlabel(ax2,'tijd'); ylabel(ax1,'m/d');
title( ax2,'Voortschrijdende waterbalans regionale aquifer');

hpos1 = area(ax1,tne(:,1),QPOS1');
hneg1 = area(ax1,tne(:,1),QNEG1');
hpos2 = area(ax2,tne(:,1),QPOS2');
hneg2 = area(ax2,tne(:,1),QNEG2');

% match the yaxes
yLim = [get(ax1,'yLim'); get(ax2,'yLim')];
yLim = [min(yLim(:,1)) max(yLim(:,2))];
set(ax1,'yLim',yLim);
set(ax2,'yLim',yLim);

iColClr = 7;
for i=1:numel(stress)
    set(hpos1(i),'faceColor',LABELS{Istress(i),iColClr},'edgeColor','none');
    set(hneg1(i),'faceColor',LABELS{Istress(i),iColClr},'edgeColor','none');
    set(hpos2(i),'faceColor',LABELS{Istress(i),iColClr},'edgeColor','none');
    set(hneg2(i),'faceColor',LABELS{Istress(i),iColClr},'edgeColor','none');
end

legend(ax1,stress{:},stress{:});
legend(ax2,stress{:},stress{:});

datetick(ax1,'x');
datetick(ax2,'x');

