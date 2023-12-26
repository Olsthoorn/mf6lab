%% mf_analyze, print and plot results fo simulation

%% Analyzing output of xsecton model
% Plot /shows the results of the smulation for a given stress period and
% time step to be set by the user
% TO 090223
clear variables
close all

%% load model name and basename contained in name.mat
load name
load(basename);  % this yields the stored matrices that make up the grid and the model
load underneath

%% load the unformatted files with the heads, the concentration and the budget terms
H=readDat([basename,'.hds']);    H=maskHC(H,[-Inf,1000],[NaN NaN]);
C=readMT3D('MT3D001.UCN');       C=maskHC(C,[ 0  ,Inf ],[NaN,NaN]);

B=readBud([basename,'.bgt'],'FLOWRIGHTFACE');
B=mf_Psi(B);

hrange = ContourRange(H,30);       dPhi = median(diff(hrange));
prange = ContourRange(B,30,'Psi'); dPsi=median(diff(prange));
crange = ContourRange(C,30);

%% Photo to plot on
foto='Series1_05_p18.jpg';
load(['.' filesep 'Photos' filesep 'photostruct.mat']); % yields struct d
i=strmatchi(foto,{d.name});
fig=ImageModel(['Photos' filesep foto],d(i).uMdl,d(i).vMdl,d(i).DX,d(i).DY); 

hold on;

% Stream function mask, to not plot stream lines above sand body
PSIMASK =  gr.size;
PSIMASK(inpolyz(gr.XM,gr.ZM,xSand,zSand)) = 1;
PSIMASK(inpolyz(gr.XM,gr.ZM,xCanL,zCanL)) = NaN;

videoName = [basename, ' March 1905'];
vidObj = VideoWriter(videoName);
vidObj.Quality = 80;
vidObj.FrameRate = 5;
vidObj.open;


for i=1:length(C)
    if i==1
        [~,hc]=contourf(gr.xc,gr.zc,XS(C(i).values),crange([1 end]),'edgeColor','none');
        [~,hh]=contour (gr.xc,gr.zc,XS(H(i).values),hrange([1 end]),'r');
        [~,hp]=contour (gr.xp,gr.zp,B(i).Psi,'w');

        xlabel('x [m]');                % set x-label
        ylabel('z [m]');                % set y-label
        set(gca,'color',[1 1 0.8])      % light blue background color
        set(gca,'clim',crange([1 end]));
    else
    
        set(hh,'zdata',XS(H(i).values));
        set(hc,'zdata',XS(C(i).values));  
        set(hp,'zdata',B(i).Psi);

        title(sprintf('%s, time = %g min',videoName,C(i).time));

        vidObj.writeVideo(getframe(gcf));
    end
end

vidObj.close();
