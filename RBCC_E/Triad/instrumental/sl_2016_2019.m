sl=cell(3,1);sl_r=cell(3,1);sl_s=cell(3,1);sl_raw=cell(3,1);
dsum=cell(3,1);
dsum_r=cell(3,1);
brewer=[157,183,185]

%read_config_
load('events_16_19.mat')



for i=1:3
    sl{i}=[];
    sl_s{i}=[];
    sl_r{i}=[];
    sl{i}.tabla=[];
    sl{i}.events={};
    sl_raw{i}=[];
    dsum_r{i}=[];
    dsum{i}=[];
for ano=2016:2019
    j=ano-2015;
   s1_=(strrep( strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Instrumental/IZO#157_sl_rw.mat','2019',num2str(ano)),'157',num2str(brewer(i))))
   if exist(s1_)
     s=load(s1_);
     sl_raw{i}=[sl_raw{i};s.sl_rw];
     sl{i}.tabla=[sl{i}.tabla;s.tabla_tc.data]
     %sl{i}.data=[sl{i}.data;s.tabla_tc.data];
     sl_s{i}=[sl_s{i},s.tabla_tc.sl];
     sl_r{i}.data=[sl_s{i},s.tabla_tc.sl_r];
     sl{i}.events=[sl{i}.events,s.tabla_tc.events];
     
     s1_=(strrep( strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/summary_old_Brw157_2019.txt','2019',num2str(ano)),'157',num2str(brewer(i))))
     s2_=(strrep( strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/summary_Brw157_2019.txt','2019',num2str(ano)),'157',num2str(brewer(i))))
 
     ds=load(s1_);
     ds_r=load(s2_);
     dsum_r{i}=[dsum_r{i};ds_r];
     dsum{i}=[dsum{i};ds];
     
    
   end
  
end
try
   t_temperature=array2table(sl{i}.tabla,'VariableNames',[{'date'},str2name((s.tabla_tc.data_lbl))],'RowNames',varname(sl{i}.events));
catch
     t_temperature=array2table(sl{i}.tabla,'VariableNames',[{'date'},str2name((s.tabla_tc.data_lbl))]);
     warning('wrong events')
end
   t_temperature.Date_str=datestr(t_temperature.date,'yyyy/mm/dd');
   t_temperature=t_temperature(:,[end,1:end-1])
   writetable(t_temperature,strrep('temp_157.xls','157',num2str(brewer(i))),'WriteRowNames',true)
   
   %%
   hold off;figure(i*4+1);for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[1,end]),'.'), end; grid; title(num2str(brewer(i)));
   axis tight; datetick('keeplimits')
   vline_v(sl{i}.tabla(:,1),'-',sl{i}.events)
   figure(i*4+2);for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[2,end]),'.'), end; grid; title(num2str(brewer(i)))
   
   figure(i*4+3);for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[1,end]),'.'), end; grid; title(num2str(brewer(i)))
    axis tight; datetick('keeplimits')
   vline_v(sl{i}.tabla(:,1),'-',sl{i}.events)
  
   figure(i*4+2);for j=1:length(sl_s{i}), hold all, ploty(sl_s{i}{j}(:,[2,end]),'o'), end; grid; title(num2str(brewer(i)))
  
   %%
   figure(i*4+4)
   plot(dsum{i}(:,1),dsum{i}(:,end),'r-',dsum{i}(:,1),dsum{i}(:,end-1),'k-')
   hold on
   plot(dsum{i}(:,1),dsum{i}(:,end),'m-',dsum_r{i}(:,1),dsum_r{i}(:,end-1),'b.')
   vline_v(events{i}.dates,'-',events{i}.labels)
   legend({'R6 op','op_ref','R6 chk','chk_ref'})
   
end

%% load configuration
%read_config_
%load alt_cfg_16_19
%load op_cfg_16_19
%events{i}=op_cfg{1}.Properties.RowNames;


%read_config_
load('events_16_19.mat')

%events{i}=op_cfg{1}.Properties.RowNames;
