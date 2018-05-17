% options_pub.outputDir=fullfile('..','html','Brewer157'); options_pub.showCode=false;
% close all; publish(fullfile(pwd,'comp_157.m'),options_pub);

%% Chequeamos eventos y corrección por SL -> creamos summarios con flag_sl=1
close all

comp.analyzed_brewer=[157 183]; % Instrumentos con los que vamos a trabajar
comp.brws_idx=find(ismember(Cal.brw,comp.analyzed_brewer));
comp.reference=find(Cal.brw==183);  % Instrumento usado como referencia

%%
for i=Cal.n_ref
    cal{i}={}; summ_orig{i}={}; summ_orig_old{i}={};
    [cal{i},summ_orig{i},summ_orig_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,...
                                          'flag_sl',1,'plot_sl',1,'flag_sl_corr',SL_corr_flag);
    % filter correction
    [summ_old{i} summ{i}]=filter_corr(summ_orig,summ_orig_old,i,A,F_corr{i});
end

Cal.n_inst=1; [ratios ratios_SL]=report_comp(Cal,comp,'summary',summ_old); 

%% Plotting
close all
 for pp=1:length(ratios)
    %%
    comp.brws_idx_=comp.brws_idx;
    figure; s={}; s_SL={};  grid;
    for ii=1:length(comp.brws_idx_)
        [aux,s_SL{comp.brws_idx(ii)}]=mean_smooth(ratios_SL{pp}(:,end),ratios_SL{pp}(:,ii+1),.125);
        [aux,s{comp.brws_idx(ii)}]=mean_smooth(ratios{pp}(:,end),ratios{pp}(:,ii+1),.125);
    end
    fprintf('\r\nOzone Deviation to %s: Operative config.\n(from %s to %s)\n',...
                           Cal.brw_name{comp.reference},datestr(ratios{pp}(1,1),2),datestr(ratios{pp}(end,1),2));   
    comp.brws_idx_(comp.brws_idx_==comp.reference)=[];
    
    h1=plot_smooth_(s{comp.brws_idx_},s_SL{comp.brws_idx_}); grid; box on;

    set(h1,'LineWidth',2); set(h1(h1~=0),'LineStyle','-');  
    title(sprintf('Relative Diffs. (%%) to the RBCC-E Triad mean\r\nDay %d%d to %d%d',...
                           diaj(min(ratios{pp}(:,1))),year(min(ratios{pp}(:,1)))-2000,...
                           diaj(max(ratios{pp}(:,1))),year(max(ratios{pp}(:,1)))-2000));     
    legend(h1, Cal.brw_str(comp.brws_idx_),strcat(Cal.brw_str{comp.brws_idx_},' SL'));
    snapnow                                           
 end
