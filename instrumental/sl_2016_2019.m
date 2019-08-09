sl=cell(3,1);sl_cr=cell(3,1);sl_raw=cell(3,1)
brewer=[157,183,185]
for i=1:3
    sl{i}=[];
    sl{i}.data=[];
    sl{i}.events={};
    sl_cr{i}=[];
for ano=2016:2019
    j=ano-2015;
   s1_=(strrep( strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Instrumental/IZO#157_sl_rw.mat','2019',num2str(ano)),'157',num2str(brewer(i))))
   if exist(s1_)
     s=load(s1_);
     sl_raw{i}=[sl_raw{i};s.sl_rw];
     %sl{i,j}=s.tabla_tc;
     sl{i}.data=[sl{i}.data;s.tabla_tc.data];
     sl{i}.events=[sl{i}.events,s.tabla_tc.events];
     
   end
  
end

   t_temperature=array2table(sl{i}.data,'VariableNames',[{'date'},str2name((s.tabla_tc.data_lbl))],'RowNames',varname(sl{i}.events));
   t_temperature.Date_str=datestr(t_temperature.date,'yyyy/mm/dd');
   t_temperature=t_temperature(:,[end,1:end-1])
   writetable(t_temperature,strrep('temp_157.xls','157',num2str(brewer(i))),'WriteRowNames',true)


end