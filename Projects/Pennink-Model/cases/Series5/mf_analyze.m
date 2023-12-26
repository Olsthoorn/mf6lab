%% mf_analyze, print and plot results fo simulation
% TO 090223 100718

clear variables
close all

load name
load(basename);  % this yields the stored matrices that make up the grid and the model
load underneath

%% load the unformatted files with the heads, the concentration and the budget terms

H=readDat([basename,'.HDS']); H=maskHC(H,[-Inf,1000],[NaN,NaN]);

C=readMT3D('MT3D001.UCN'   ); C=maskHC(C,[-90,Inf],[NaN,NaN]);

B=readBud([basename,'.bgt'],'FLOWRIGHTFACE');
B=mf_Psi(B);  for it=1:numel(B), B(it).Phi(IBOUND==0)=NaN; end

hrange = ContourRange(H,30);
crange = ContourRange(C,30);
prange = ContourRange(B,30,'Psi');

%%
foto='Series5_01_p70.jpg';
load(['.' filesep 'Photos' filesep 'photostruct.mat']); % yields struct d
it=strmatchi(foto,{d.name});
fig=ImageModel(['Photos' filesep foto],d(it).uMdl,d(it).vMdl,d(it).DX,d(it).DY); 

lightBlue = [1 1 0.8];
set(gca','nextplot','add','color',lightBlue,'clim',crange([1 end]));
xlabel('x [m]');                 % set x-label
ylabel('z [m]');                 % set y-label

videoName=[basename, 'March 1905'];
vidObj=VideoWriter(videoName);
vidObj.FrameRate = 5;
vidObj.Quality = 80;
vidObj.open();

time = [C.time];

tts = sprintf('%s, time = %%g min',videoName);

for it=1:length(time)
    
    tts1= sprintf(tts,time(it));
    
    if it==1
           ht = title(tts1);
        [~,hc]=contourf(gr.xc,gr.zc,XS(C(it).values),crange,'edgeColor','none');
        [~,hh]=contour (gr.xc,gr.zc,XS(H(it).values),hrange,'r');
        [~,hp]=contour (gr.xp,gr.zp,B(it).Psi,prange,'w');
    else
        set(ht,'string',tts1);
        set(hc,'zdata',XS(C(it).values));
        set(hh,'zdata',XS(H(it).values));
        set(hp,'zdata',B(it).Psi);

        vidObj.writeVideo(getframe(fig));
    end
end

vidObj.close();
