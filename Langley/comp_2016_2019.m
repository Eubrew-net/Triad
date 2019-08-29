
%path_root='D:\CODE\iberonesia\RBCC_E\Triad\2019\Triad'
%path_root='D:\CODE\iberonesia\RBCC_E\2019\Triad\langley'
%if ismac
% path_root='/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/Triad/2019/Triad'
%end

run('../instrumental/read_config_')

% remove year from cal.path_root
path_root=fullfile(Cal.path_root,'..');
fblacklist=fullfile(path_root,'configs',['blacklist_',num2str(brewer(i)),'.txt']);  
confname=["Op","Alt"];
Width=20; Height=10; 

ds_o=cell(3,1);
ds_a=ds_o;
ds=cell(3,2);
brewer=[157,183,185];
ano0=2014;



for i=1:3
    ds_o{i}=[];
    ds_a{i}=[];
    for ano=2015:2019
        j=ano-2014;
        %s1_=(strrep( strrep(fullfile(path_root,'summary_Brw157_2019.txt'),...
        %    '2019',num2str(ano)),'157',num2str(brewer(i))))
        s1_=fullfile(path_root,'Triad',num2str(ano),'Triad',['summary_Brw',num2str(brewer(i)),'_',num2str(ano),'.txt']);
        %s2_=(strrep( strrep(fullfile(path_root,'summary_old_Brw157_2019.txt'),...
        %    '2019',num2str(ano)),'157',num2str(brewer(i))));
        s2_=fullfile(path_root,'Triad',num2str(ano),'Triad',['summary_old_Brw',num2str(brewer(i)),'_',num2str(ano),'.txt']);
       
        if exist(s1_)
            s=load(s1_);
            ds_a{i}=[ds_a{i};s(:,1:15)];  % write summary incorpora ahora A1 y ETC para recalcular el ozono, los ficheros antiguos no lo hacen
        else
           warning([s1_ ,'not found'])
        end    
        
        if exist(s2_)
            s=load(s2_);
            ds_o{i}=[ds_o{i};s(:,1:15)];
        else
         warning([s2_ ,'not found'])
        end
        
        fprintf('Año %d brewer %d diff(o1,a2)=%d diff(o1,a2)=%d\n',ano,brewer(i),sum((ds_o{i}(:,6)-ds_a{i}(:,10))~=0),sum((ds_o{i}(:,10)-ds_a{i}(:,6))~=0));
    end
    
    %%Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc
    
 
    
 ds{i,1}=blacklist_summary(fblacklist,ds_o{i});
 ds{i,2}=blacklist_summary(fblacklist,ds_a{i});
     
end

%%
for j=1:2
 
TSYNC=10; % (min)
ref1=fix(ds{1,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
ref2=fix(ds{2,j}(:,1)*24*60/TSYNC)/24/60*TSYNC; 
ref3=fix(ds{3,j}(:,1)*24*60/TSYNC)/24/60*TSYNC;
s=scan_join([ref1,ds{1,j}],[ref2,ds{2,j}]);
s=scan_join(s,[ref3,ds{3,j}]);
label{1}='Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc';
% s 1 date
%        2:16  157 'Date sza m2 temp nd O3_1 std ms9_corr ms9 O3_2 std O3_1_sl std R6_ref R6_calc'
%       17:31  183
%       32:45  185

% datos simultaneos
jsim=all(~isnan(s(:,7:15:end)),2) ; 
ref=median(s(jsim,7:15:end),2);
sim=s(jsim,:);

% prueba
%mydata=sim(:,7:15:end)-mean(sim(:,7:15:end),2);
%figure
%plot(sim(:,1),mydata)
%prueba

% medidas que difieren en menos de 1º
%histogram(sim(:,3:15:end)-mean(sim(:,3:15:end),2))
jsim2=all(abs(sim(:,3:15:end)-mean(sim(:,3:15:end),2))<1.0,2);

ref1=median(sim(jsim2,7:15:end),2);  % reference 1st configuration  (O3_1)
ref2=median(sim(jsim2,11:15:end),2);  % reference 2st configuration  (O3_2)
ref3=median(sim(jsim2,13:15:end),2);  % reference 1st configuration+SL  (O3_1_sl)

ratio1=[sim(jsim2,1),100*(sim(jsim2,7:15:end)-ref1)./ref1];    % 1st configuration
ratio2=[sim(jsim2,1),100*(sim(jsim2,11:15:end)-ref2)./ref2];   % 2nd  configurationn
ratio3=[sim(jsim2,1),100*(sim(jsim2,13:15:end)-ref3)./ref3];   % 1st configuraiton + SL correction 

m=[sim(jsim2,1),sim(jsim2,4:15:end)]; % airmas

x1=smoothdata(ratio1(:,2:end),'gaussian',15,'SamplePoints',ratio1(:,1));
x2=smoothdata(ratio2(:,2:end),'gaussian',15,'SamplePoints',ratio2(:,1));
x3=smoothdata(ratio3(:,2:end),'gaussian',15,'SamplePoints',ratio3(:,1));

f=figure;
set(f,'Tag','comp_2016_2019');
plot(ratio1(:,1),x1,'.'); hold all
%plot(ratio2(:,1),x2,'-','LineWidth',2); hold all
%plot(ratio2(:,1),x3,':'); hold all
grid;
datetick('x',12);
%title(['Configuracion',num2str(j)])
title(strcat('Conf. ',confname(j)))
set(gca,'Ylim',[-2,2]);
end

ix=sort(findobj('-regexp','Tag','comp_2016_2019'));

printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);
