function [xx,yy]=MeasureOnPhoto(str)
% [xx,yy[=MeasureOnPhoto(str)
% Select anything on phote from list in photostruct.mat using searchstr str
% by clicking. ENTER to stop. The coordinates are then given.
%
% Olsthoorn 100721 100727

load('photostruct');

if nargin==0,
    fprintf('Select one of the following photos using unique combination of first letters\n');
    fprintf('for instance:\n');
    fprintf('[xx,y]=MeasureOnPhoto(''Series_01'')\n');
    fprintf('When the crossing hairlines appear on the photo select\n');
    fprintf('your coordinates. Stop by pressing enter.\n');
    fprintf('The selected coordinates will be in the output of the function.\n');
    fprintf('Available photos are:\n');
    for i=1:length(d)
       fprintf('PhoteNr %2d, %-20s u=[%3d %3d]  v=[%3d %3d]  Dx=%.1f Dz=%.1f cm\n',...
           i,d(i).name,round(d(i).uMdl),round(d(i).vMdl),d(i).DX,d(i).DY);
    end
    fprintf('Now select your points for photo Series_01_p04.jpg:\n');
    str='Series1_01';
end

i=strmatchi(str,{d.name});

fig=ImageModel(d(i).name,d(i).uMdl,d(i).vMdl,d(i).DX,d(i).DY); hold on;
xlabel('x cm'); ylabel('z cm'); title(sprintf('photo %s, click [ENTER=stop]',d(i).name));
[xx,yy]=ginput;
