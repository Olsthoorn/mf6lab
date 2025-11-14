
for i=1:Nsec;
GHGn_mv = GHGn-AHN;% reken GHG numeriek van NAP naar m tov mv
GVGn_mv = GVGn-AHN; 
GLGn_mv = GLGn-AHN; 

GHGa_mv = GHGa-AHN;% reken GHG analytisch van NAP naar m tov mv
GVGa_mv = GVGa-AHN; 
GLGa_mv = GLGa-AHN; 
 
% P(i).GHG = ghg(i);
% P(i).GVG = gvg(i); 
% P(i).GLG = glg(i);
end
% [P.GHG]=[P.GHG]-[P.z0];
% [P.GVG]=[P.GVG]-[P.z0];
% [P.GLG]=[P.GLG]-[P.z0];

% %%%tussen=repmat([P.FID2]',1,size(h,2));
% name= reshape(repmat([P.FID2]',1,size(h,2))',1,[])';
% filternr=repmat(ones(size(h,2),1),size(h,1),1);
% input_menyanthes=[name,filternr,(inp_menyanthes(:,1)-repmat(datenum(1899,12,30),1,size(inp_menyanthes(1)))),inp_menyanthes(:,2)]; %tbv van filternr in data sheet
% %%%filterno=ones(size(h,1),1);% tbv filternr in gegevens sheet
% time= reshape(repmat(t,1,size(h,1)),1,[])';


xlswrite('ggor_Noorderpark_2016_kw3.xlsx',[[FID],[b],[D1],[AHN],[c],[GLGa_mv],[GVGa_mv],[GHGa_mv],[GLGn_mv],[GVGn_mv],[GHGn_mv],[b],[c],[vk1],[hk1],[sy1],[GP]],'results_GGOR','A2');
xlswrite('ggor_Noorderpark_2016_kw3.xlsx',[[bodcod1]],'results_GGOR','R2');
%xlswrite('HELP2005toets agor_Noorderpark',[[P.GHGanalytic]',[P.GLGanalytic]'],'bodemsoort','C2');% gebruik de analytische oplossing

%% schrijf data naar de HELP tabellen
% xlswrite('HELP2005toets_agor_Noorderpark_flex2015_kw1.xls',[[P.GHG]',[P.GLG]'],'bodemsoort','C2');% gebruik de numerieke oplossing
% for i=1:length(P),P(i).BODCOD1=deblank(P(i).BODCOD1); end 
% xlswrite('HELP2005toets_agor_Noorderpark_flex2015_kw1.xls',[[P.BODCOD1]'],'bodemsoort','B2');
% xlswrite('HELP2005toets_agor_Noorderpark_flex2015_kw1.xls',[[P.FID3]'],'bodemsoort','A2');
