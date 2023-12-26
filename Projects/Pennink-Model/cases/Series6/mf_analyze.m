%% mf_analyze, print and plot results fo simulation
% Pennink Series6
%
% TO 090223 100721 100726 120517

clear variables
close all

%% load model name and basename contained in name.mat
load name
load(basename);  % this yields the stored matrices that make up the grid and the model
load underneath
mf_checkdir;

%% load the unformatted files with the heads, the concentration and the budget terms

H=readDat([basename,'.HDS']); H=maskHC(H,[-Inf,1000],[NaN,NaN]);

C=readMT3D('MT3D001.UCN');    C=maskHC(C,[-90,Inf ],[NaN,NaN]);

B=readBud([basename,'.BGT']);
B=mf_Psi(B);

crange=ContourRange(C,50);
hrange=ContourRange(H,50);
prange=ContourRange(B,50,'Psi');

% Stream function mask, to not plot stream lines above sand body
PSIMASK = ~(inpolyz(gr.XM, gr.ZM ,xSand,zSand) | ...
            inpolyz(gr.XM, gr.ZM ,xCanL,zCanL) | ...
            inpolyz(gr.XM, gr.ZM ,xCanR,zCanR));

for i=1:numel(B), B(i).Psi(PSIMASK)=NaN; end

%% Set up contours and first plot to initialize animation

fig=figure('position',screenPos(0.80));

ax = axes('nextplot','add','cLim',crange([1 end]),'parent',fig);

sandcolor = [1 1 0.6];
xlabel(ax,'x [m]');
ylabel(ax,'z [m]');

ts1 = sprintf('Pennink (1915) sand model, experiment series 6 (1905); time = %%g min');

time = [C.time];

vidObj=VideoWriter(basename);
vidObj.FrameRate = 3;
vidObj.Quality = 80;
vidObj.open;

for it=1:length(time)
    
    ts2 = sprintf(ts1,time(it));

    if it==1
        patch(xSand,zSand,sandcolor,'parent',ax);

        ht = title(ts2,'fontsize',15);
        [~,hc]= contourf(ax,gr.xc,gr.zc,XS(C(it).values),crange,'edgecolor','none');
        [~,hh]= contour( ax,gr.xc,gr.zc,XS(H(it).values),hrange,'r');
        [~,hp]= contour( ax,gr.xp,gr.zp,B(it).Psi,prange,'edgeColor',grey);
        
        patch(xUnsat,zUnsat,sandcolor,'parent',ax);

        patch(xCanL,zCanL,'w','parent',ax);
        patch(xCanR,zCanR,'w','parent',ax);

        %patch(xCanM,zCanM,'w');
        patch(xWatL,zWatL,'b','parent',ax);
        patch(xWatR,zWatR,'b','parent',ax);
        
        patch(xWell,zWell,'w','parent',ax);

        mf_logo(ax);
        
        axis('tight');
    else
        set(ht,'string',ts2);
        set(hc,'Zdata',XS(C(it).values));
        set(hh,'Zdata',XS(H(it).values));
        set(hp,'Zdata',B(it).Psi);
    end
        
    vidObj.writeVideo(getframe(gcf));
end

vidObj.close;

%% Discharge through boundary conditions

fprintf('%10s %10s %10s %10s %10s %10s %10s\n','Period','tstep','QL','QM','QR','QS','Qtot');
ccmpm2lph=1e-6*60; % m3/min to L / h
for it=1:length(B)
    iGHB = strmatchi('HEADDEPBOUNDS',B(it).label);
    iCHD = strmatchi('CONSTANTHEAD' ,B(it).label);

    if iGHB
        QL=sum(B(it).term{iGHB}(IBOUND==iWatL))*ccmpm2lph; % in liters per hours
    else
        QL=0;
    end
    if iCHD
        QM=sum(B(it).term{iCHD}(IBOUND==iWatR))*ccmpm2lph;
        QR=sum(B(it).term{iCHD}(IBOUND==iWell))*ccmpm2lph;
        QS=sum(B(it).term{iCHD}(IBOUND==iMilk))*ccmpm2lph;
    else
        QM = 0;
        QW = 0;
        QS = 0;
    end
    fprintf('%10d %10d %10g %10g %10g %10g %10g\n',B(it).period,B(it).tstp,QL,QM,QR,QS,sum([QL QM QR QS]));
end

%% Interface computation
%  The head in the milk is set in the input. It serves to keep the
%  interface at the desired elevation. However, the correct milk-head is
%  known only afterwards and iteration is cumbersum. Therefore, we compute the
%  desired milk-head given the freshwater head and the saltwater head at
%  two strategic positions in the model (below the right ditch at z=50 cm,
%  which should be a good position for the freshwater head without much
%  influence of the converging flow near te ditch and still fresh, and the
%  point near the milk supply. If the correct head is chosen, the interface
%  should be stable and little milk should flow and so this saltwater
%  head should be more or less that in the milk reservoir).
rhof=1000; rhos=1029; zIface=23;  % must match what has been specified in model

I=xyzindex([95 95]',[0 0]',[50 5]',gr); % point positions

phif=NaN(size(H));  % freshwater head in first  point
phis=NaN(size(H));  % saltwater  head in second point 
for it=1:length(H)
    phif(it)=H(it).values(I(1));
    phis(it)=H(it).values(I(2));
end
figure; hold on

subplot(2,1,1);  % plot phif, phis and phis desired to maintain zIface
plot(time,phif,'b',...
     time,phis,'r',...
     time,(rhos-rhof)/rhos*zIface+rhof/rhos.*phif,'g');
legend('phi_f(50)','phi_s(5)','phi_s desired');
xlabel time; ylabel head; title('fresh and saltwater head and satwater head need to keep interface at zIface');

subplot(2,1,2); % plot computed interface form heads and desired zIface
plot(time,(rhos*phis-rhof*phif)/(rhos-rhof),'r',...
     time,zIface*ones(1,length(H)),'g');
legend('computed interface elevation','desired elevation zIface')
xlabel time; ylabel head; title('computed interface given the heads in the model');
