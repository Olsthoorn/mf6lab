%% mf_analyze, print and plot results fo simulation

%% Analyzing output of xsecton model
% Plot /shows the results of the smulation for a given stress period and
% time step to be set by the user
% TO 090223 111212

clear variables
close all

film=1;  % if film==1, make avi file else not

%% load model name and basename contained in name.mat
load name
load(basename);  % this yields the stored matrices that make up the grid and the model
load underneath

%% Start

%% load the unformatted files with the heads, the concentration and the budget terms

H=readDat([basename,'.hds']); H = maskHC(H,[-Inf,1000],[NaN,NaN]);

C=readMT3D('MT3D001.UCN');    C = maskHC(C,[-90,1000],[NaN,NaN]);

B=readBud([basename,'.bgt'],'FLOWRIGHTFACE');
B=mf_Psi(B);

prange=ContourRange(B,30,'Psi');   dPsi = median(diff(prange));
hrange=ContourRange(H,30);         dPhi = median(diff(hrange));
crange=ContourRange(C,50);

% Stream function mask, to not plot stream lines above sand body
% PSIMASK=NaN(Nz+1,Nx-1);
% PSIMASK(inpolygon(ones(size(zGr(:)))*xGr(2:end-1),zGr(:)*ones(size(xGr(2:end-1))),xSand,zSand))=1;
% PSIMASK(inpolygon(ones(size(zGr(:)))*xGr(2:end-1),zGr(:)*ones(size(xGr(2:end-1))),xCanL,zCanL))=NaN;

%% Photo to plot on
foto='Series3_02_p52.jpg';
load(['.' filesep 'Photos' filesep 'photostruct.mat']); % yields struct d
i=strmatchi(foto,{d.name});
fig=ImageModel(['Photos' filesep foto],d(i).uMdl,d(i).vMdl,d(i).DX,d(i).DY); 

%set(fig,'position',get(0,'ScreenSize'));
lightBlue = [1 1 0.8];

set(gca,'nextplot','add','color',lightBlue,'clim',crange([1 end]));
xlabel('x [m]');                 % set x-label
ylabel('z [m]');                 % set y-label

videoName = [basename, 'March 1905'];
vidObj    = VideoWriter(videoName);
vidObj.FrameRate = 5;
vidObj.Quality   = 80;
vidObj.open();

for i=1:length(C);
    if i==1
        [~,hc]=contourf(gr.xc,gr.zc,XS(C(1).values),crange([1 end]),'edgeColor','none');
        [~,hh]=contour (gr.xc,gr.zc,XS(H(1).values),hrange([1 end]),'r');
        [~,hp]=contour (gr.xp,gr.zp,B(i).Psi,prange([1 end]),'w');

    else    
        set(hh,'zdata',XS(H(i).values));
        set(hc,'zdata',XS(C(i).values));   
        set(hp,'zdata',B(i).Psi);
    end
    
    set(get(hc,'children'),'edgecolor','none');
    
    title(sprintf('%s, time = %g min',videoName,C(i).time));

    vidObj.writeVideo(getframe(gcf));
end
