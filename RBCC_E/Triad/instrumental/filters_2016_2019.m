%% filter analysis
% this file is CODE/iberonesia/RBCC_E/Triad/instrumental/filters_2016_2019.m
%
% auxiliary functions included at the end of this file: 
%    handles2save = changeFigures(Cal, figure_tags, fig_xlim, events)
%    save_figs(Cal, figure_handles, file_suffix, fig_width, fig_height)
%    [summary, summary_old] = load_summ(brw_idx, Cal)
%    [config, ozone_ds, ozone_raw] = load_bdata(Cal)
%    hline_pos
%
% 2DO:
%   
%
% 201909 JLS, based on code from everybody


%% input options

% ----------------------------
% start and end analysis dates
date_start=datenum(2016,1,1);
date_end=now;

% ----------------------------
brewers2check=[1,2,3]; % 1->157, 2->183, 3->185 

% ----------------------------
verbose_output=false; % print some more data on the command window

% ----------------------------
% get the ETC filter corrections from the FIOAVG file
filters_fioavg_run=true; % run this analysis? true/false

% for filters_fioavg, we will save the figure "FI_TIME_ETC2" from the whole set of figures generated:
%     in the call to filter_rep:
%         Attenuation Filter Test, Difference with respect to nominal values -> tag="FI_STATS"
%         Attenuation Filter Test, Diff. (%) with respect to nominal values [...] -> tag="FIOSTATS_1" (nominal), "FIOSTATS_2" (operative)   
%         Wavelength dependence of the attenuation filter -> tag="FI_wavelength"
%     in the call to filters_data:
%         ND attenuations vs. time (sample slits #3, top, and #5, bottom) -> tag="FI_TIME_atts"
%         Filter ETC correction vs time. Monthly means -> tag="FI_TIME_ETC2"

% ----------------------------
% get the ETC filter corrections from the ozone measurements
filters_ozone_run=true; % run this analysis?

% for filters_ozone, we will save the figures "FILTER_DISTRIBUTION", and "Ozone_diff_filter_rel_comparative" from the whole set:
%     in the call to ozone_filter_analysis_mi:
%         percentage of measurements with each filter (pie chart) -> tag="FILTER_DISTRIBUTION"
%         ozone difference by filter change -> tag="Ozone_diff_filter"
%         ozone relative difference by filter change -> tag="Ozone_diff_filter_rel"
%     after the call to ozone_filter_analysis_mi:
%         ozone relative difference by filter change, summ_old vs summ -> tag="Ozone_diff_filter_rel_comparative"
ozone_filter_analysis_mi_figs2save=[];

% ----------------------------
% get the ETC filter corrections looking at the langleys
filters_langley_run=true; % run this analysis?

% for the filters_langley analysis, load the config, ozone_ds, and ozone_raw data from the file in filters_langley_data_file?
% if false, the data will be generated from the B files and saved to filters_langley_data_file
filters_langley_load_data=true;
filters_langley_data_file=["filters_2016_2019.mat"];

% for the filters_langley analysis, do an AOD screening using an AERONET file (such as eg CODE/iberonesia/RBCC_E/2019/Triad/BSRN/190101_191231_Izana.lev15)
% leave empty ([]) for no screening
filters_langley_aod_file=[];

% for the filters_langley analysis, do a cloud screening using a file with data from the BSRN (such as eg CODE/iberonesia/RBCC_E/2019/Triad/BSRN/cloudScreening.txt)
% leave empty ([]) for no screening
filters_langley_cloud_file=[];

% for the filters_langley analysis, the following figures will be saved:
%     PDF vs ETC, for all filters, per event -> tags are ETC_FILTER_HISTOGRAM_event_$event_num
%     Filter ETC correction vs time, for all filters -> tag is ETC_FILTER_SERIES

% ----------------------------
% suffix to append to the filenames of the figures and latex tables
file_suffix=char("jlsTest"); % NOTE this MUST be a char vector! so either use '...' or char("...")

% ----------------------------
% general fig options
fig_width=20;
fig_height=10;
fig_xlim=[date_start,date_end];
fig_linewidth=1;


%% define filter events for each brewer
% the order in the variable Cal.brw defines the index of each brewer
filter_events{1}.dates=[datenum(2016,1,1)]; 
filter_events{1}.labels={'Original filters'};

filter_events{2}.dates=[datenum(2016,1,1), datenum(2018,3,1)]; 
filter_events{2}.labels={'Original filters', 'Power supply change'};

filter_events{3}.dates=[datenum(2016,1,1),datenum(2017,3,1),datenum(2018,2,4)]; 
filter_events{3}.labels={'Original filters','Filter change #1','Filter change #2'};


%% **********************************************
%  **********************************************
%  **********************************************
%
%                code starts here
%
%  hopefully, you won't need look past this line!
%
%  **********************************************
%  **********************************************
%  **********************************************


%% read configs and set paths
% "read_config_" loads "calizo_setup" and defines a lot of vars, including
%    Cal.path_root=fullfile(Cal.path_root,'RBCC_E',num2str(Date.cal_year));
%    Cal.brw=[157,183,185];
run(fullfile('..','read_config_'));


%% change date range from the one defined in read_config_ 
Cal.Date.day0=date_start;
Cal.Date.dayend=date_end;
Date.CALC_DAYS=Cal.Date.day0:Cal.Date.dayend; 
Cal.Date=Date;


%% ***************************************************
%  ***************************************************
%
%  get the ETC filter corrections from the FIOAVG file
%
%  ***************************************************
%  ***************************************************

if filters_fioavg_run
    
    % loop over all the brewers
    for brw_idx=brewers2check % usally, Cal.n_ref=[1 2 3]
        close all;
        
        % initialize table
        tabla_fioavg_refF2=[];
        
        % the brewer to analyze is Cal.brw{Cal.n_inst}
        Cal.n_inst=brw_idx;
        fprintf(strcat("\n\n******** running filters_fioavg for Brewer #",Cal.brw_str{Cal.n_inst},"\n"))
              
        % run report_filter
        
        % to get nicer plots, give report_filter some more data at the beginning of the range
        %filter_events_mod=filter_events{brw_idx};
        %filter_events_mod.dates(1)=filter_events_mod.dates(1)-60;       
        %tabla_fi=report_filter(Cal,'grp_custom',filter_events_mod,'date_range',[filter_events_mod.dates(1),now]);
                
        tabla_fioavg=report_filter(Cal,'grp_custom',filter_events{brw_idx},'date_range',[filter_events{brw_idx}.dates(1),now]);
        
        % use the aux function changeFigures (see the end of this file) to add the brewer id to the figs' tags,
        %+set the x axis, and add the events
        handles2save=changeFigures(Cal, "FI_TIME_ETC2", fig_xlim, filter_events{brw_idx});
        
        % save the figures using the aux function saveFigs (see the end of this file) to call printfiles_report with
        %+appropriate parameters (in particular, an aux_pattern cell with the file_suffix defined in the input options)
        saveFigs(Cal, handles2save, file_suffix, fig_width, fig_height, fig_linewidth);
                    
        % show table with the Filter ETC corr per event
        fprintf("** ETC Filter corr, all filters by separate\n")
        displaytable(tabla_fioavg.data(:,2:2:(end-1)), tabla_fioavg.data_lbl(1:2:(end-1)), 10, '.2f', tabla_fioavg.events);
        
        % table taking as reference filter 2
        tabla_fioavg_refF2.data(:,1)=tabla_fioavg.data(:,1); % date
        tabla_fioavg_refF2.data(:,2:6)=tabla_fioavg.data(:,2:2:(end-1)); % filters will be in cols 2-6
        tabla_fioavg_refF2.data_lbl=tabla_fioavg.data_lbl(1:2:(end-1));
        tabla_fioavg_refF2.events=tabla_fioavg.events;
        tabla_fioavg_refF2.data(:,2:end)=tabla_fioavg_refF2.data(:,2:end)-tabla_fioavg_refF2.data(:,3); % filter 2 is in col 3
        
        fprintf("**** ETC Filter corr from FIOAVG, taking as reference the correction to Filter #2\n")
        displaytable(tabla_fioavg_refF2.data(:,2:end), tabla_fioavg_refF2.data_lbl(1:end), 10, '.2f', tabla_fioavg_refF2.events);
        
        % save table as a latex file
        fprintf(strcat("** saving table to ", Cal.dir_tables, "\n"))
        
        % make sure file_suffix is a vector of chars or the following line will break!
        table_file=fullfile(Cal.dir_tables, ['table_ETC_FILTER_CORR_FIOVAG_', Cal.brw_str{Cal.n_inst}, '_', file_suffix, '.tex'] );
        
        % to change rows to cols, add ' to tabla_fi.data(:,2:end), and change columnlabels to rowlabels and viceversa
        matrix2latex_ctable(tabla_fioavg_refF2.data(:,2:end), table_file, ...
                            'columnlabels', str2latex(tabla_fioavg_refF2.data_lbl), ...
                            'rowlabels', str2latex(tabla_fioavg_refF2.events), ...
                            'alignment', 'c');
                        
    end % end loop over brewers
     
end % end filters_fioavg


%% **********************************************************
%  **********************************************************
%
%  get the ETC filter corrections from the ozone measurements
%
%  **********************************************************
%  **********************************************************

if filters_ozone_run
    
    % initialize storage
    % i'm not sure what's the difference between these two summaries... do I have to check both of them?
    summ=cell(length(Cal.n_ref),1);
    summ_old=cell(length(Cal.n_ref),1);
    
    %% loop over all the brewers
    for brw_idx=brewers2check %Cal.n_ref
        close all;
        
        % initialize output table
        tabla_ozone.data=[];
        tabla_ozone.data_lbl={'F#1 corr', 'F#2 corr', 'F#3 corr', 'F#4 corr', 'F#5 corr'};
        tabla_ozone.events=filter_events{brw_idx}.labels;
        
        % the brewer to analyze is Cal.brw{Cal.n_inst}
        Cal.n_inst=brw_idx;
        fprintf(strcat("\n\n******** running filters_ozone for Brewer #",Cal.brw_str{Cal.n_inst},"\n"))
   
        %% load and format the data using the aux function load_summ (see the end of this file)
        [summ{brw_idx}, summ_old{brw_idx}]=load_summ(brw_idx, Cal);
        
        %% do the analysis per event
        num_events=length(filter_events{brw_idx}.dates);
        for event_idx=1:num_events
            fprintf(strcat("**** working on event #",num2str(event_idx),"\n"))
            
            %% get only the data between this event and the next (or the end of the file, if the last event)
            summ_event{brw_idx}=sliceVarByEvent(summ{brw_idx}, filter_events{brw_idx}.dates, event_idx, num_events);
            summ_old_event{brw_idx}=sliceVarByEvent(summ_old{brw_idx}, filter_events{brw_idx}.dates, event_idx, num_events);
            
            
            %% run ozone_filter_analysis_mi
            % a max time of 15 minutes between measurements with different filters is the default hardcoded in ozone_filter_analysis_mi
            % it might be too high, and perhaps it could be better to switch the time threshold to an airmass one, see below
            fprintf("** No. of filter changes (max time diff between measurements: 15 minutes)")
            rel1=ozone_filter_analysis_mi(summ_event,Cal,Cal.n_inst);
            %rel2=ozone_filter_analysis_mi(summ_old_event,Cal,Cal.n_inst,0);

            
            %% change the tag and title of the FILTER_DISTRIBUTION fig just generated to add the event number and date period
            h=findobj('Tag','FILTER_DISTRIBUTION');
            set(h,'Tag',strcat('FILTER_DISTRIBUTION_event_',num2str(event_idx)));
            figure(h);
            title(sprintf('%s filter utilization %% (%s to %s)',Cal.brw_name{Cal.n_inst}, ...
                datestr(summ_event{brw_idx}(1,1),29), datestr(summ_event{brw_idx}(end,1),29) ));
            
            hold off
            %% plot relative ozone differences when filters are changed 
            % (code stolen from CODE/iberonesia/RBCC_E/2019/Triad/Langley/ozone/Langley157.m)
            f=figure;  
            set(f,'Tag',strcat('Ozone_diff_filter_rel_comparative_event_',num2str(event_idx)));
            
            p=patch([0.5 8.5 8.5 0.5],[-0.5 -0.5 0.5 0.5],[.93 .93 .93]);
        
            h1=boxplotCsub(rel1(:,2:2:end),1,'o',1,1,'g',true,1,true,[1 1],1.5,0.005,false);
            %h2=boxplotCsub(rel2(:,2:2:end),1,'*',1,1,'r',true,1,true,[1 1],1.25,0.05,false);
        
            %set(gca,'YLim',[-2 2],'XTickLabel',{'0-64','64-0','64-128','128-64','128-192','192-128','192-256','256-192'});
            set(gca,'YLim',[-2 2],'XTickLabel',{'0-1','1-0','1-2','2-1','2-3','3-2','3-4','4-3'});
        
            ylabel('Relative Difference (%)');
            xlabel('Filter chg.');
            hline(0,'-.k');
            grid;
        
            title(sprintf('Ozone Relative differences by filter chg. %s\r\n Referenced always to lower filter for each group',Cal.brw_name{Cal.n_inst}));
        
            %legend([h1(end,1),h2(end,1)],{'summ','summ\_old'},'Location','SouthWest','Orientation','Horizontal');
            legend([h1(end,1)],{'summ'},'Location','SouthWest','Orientation','Horizontal');
          
            %% ETC filter corr from ozone changes
            % I guess only summ (and not both summ and summ_old) have to be checked?
            summ_event_mod=summ_event;
            
            % Delta(ETC filter corr) = -mu*alpha*Delta(ozone) if the ozone has been calculated without the filter
            % correction (otherwise, Delta(ozone) should be zero!)
            % So if we call ozone_filter_analysis_mi with the ozone multiplied by the airmass and the abs coeff, we will 
            % get the Delta(ETC filter corr)
            %summ_event_mod{brw_idx}(:,6)=summ_event_mod{brw_idx}(:,6).*summ_event_mod{brw_idx}(:,3).*summ_event_mod{brw_idx}(:,16)*10;
            %[~, abs1]=ozone_filter_analysis_mi(summ_event_mod, Cal, Cal.n_inst, 0);     

            % However, it is also true that Delta(ETC filter corr) = -Delta(MS9), and this seems safer because the MS9
            % are always free of filter corrections (I think...). In this case, instead of ozone_filter_analysis_mi, use 
            % ozone_filter_analysis_mi_ms9 to get differences of ms9 (instead of differences of ozone) 
            %
            % the _ms9 function also introduces another difference: instead of removing data from filter changes that took
            % place far apart in time, it removes them if the diff in airmass is larger than some value. in this case, it seems
            % more natural to use the airmass than the time, because the airmass enters directly in the formula, and indeed
            % it has been easier to get reasonable results playing with an airmass threshold than with a time one.
            % to analyze the diffs in time and airmass when the filters are changed, you can use e.g.
            % [datestr(summ_event_mod{3}(:,1)), string(summ_event_mod{3}(:,3)), string(summ_event_mod{3}(:,5))]
            fprintf("** No. of filter changes (max airmass diff between measurements: 0.03)")
            [~, abs1]=ozone_filter_analysis_mi_ms9(summ_event_mod, Cal, Cal.n_inst, 0);
                        
        
            % filter 0->1 delta
            % ozone_filter_analysis_mi splits the data in 0->1 and 1->0, so join it
            % differences from ozone_filter_analysis_mi are always wrt the lowest filter
            filter_0_1=cat(1, abs1(:,2), abs1(:,4) );
            delta_etc_filter_0_1=-nanmedian(filter_0_1);
        
            filter_1_2=cat(1, abs1(:,6), abs1(:,8) ); 
            delta_etc_filter_1_2=-nanmedian(filter_1_2);
        
            filter_2_3=cat(1, abs1(:,10), abs1(:,12) ); 
            delta_etc_filter_2_3=-nanmedian(filter_2_3);
        
            filter_3_4=cat(1, abs1(:,14), abs1(:,16) ); 
            delta_etc_filter_3_4=-nanmedian(filter_3_4);
        
            delta_etc_filter_2_4=delta_etc_filter_2_3+delta_etc_filter_3_4;
            
            % output the calculated values now instead of in table later
            if verbose_output
                fprintf(strcat("** ETC Filter corrections for event ", num2str(event_idx), "\n"))
                fprintf(strcat("0->1: ",num2str(delta_etc_filter_0_1),"\n"))
                fprintf(strcat("1->2: ",num2str(delta_etc_filter_1_2),"\n"))
                fprintf(strcat("2->3: ",num2str(delta_etc_filter_2_3),"\n"))
                fprintf(strcat("3->4: ",num2str(delta_etc_filter_3_4),...
                    " 2->4: ", num2str(delta_etc_filter_2_4),"\n"))
            end
            
            % save data of this event to a table, to save it to a file at the end of the events' loop
            tabla_ozone.data(event_idx,1)=filter_events{brw_idx}.dates(event_idx);
            tabla_ozone.data(event_idx,2)=-delta_etc_filter_1_2; % filter #1 wrt filter #2
            tabla_ozone.data(event_idx,3)=0; % filter #2 is the reference
            tabla_ozone.data(event_idx,4)=delta_etc_filter_2_3; % filter #3
            tabla_ozone.data(event_idx,5)=delta_etc_filter_2_4; % nothing for filter #4
            tabla_ozone.data(event_idx,6)=NaN; % nothing for filter #5 because it is not even considered in the calc!
            
        end % loop over events
        
        % show again the Filter ETC corr per event but this time as a table
        fprintf("**** ETC Filter corr from ozone measurements, taking as reference the correction to Filter #2\n")
        displaytable(tabla_ozone.data(:,2:end), tabla_ozone.data_lbl, 10, '.2f', tabla_ozone.events);
        
        % save table
        fprintf(strcat("** saving table to ", Cal.dir_tables, "\n"))
        
        % make sure file_suffix is a vector of chars or the following line will break!
        table_file=fullfile(Cal.dir_tables, ['table_ETC_FILTER_CORR_OZONE_', Cal.brw_str{Cal.n_inst}, '_', file_suffix, '.tex'] );
        
        matrix2latex_ctable(tabla_ozone.data(:,2:end), table_file, ...
                        'columnlabels', str2latex(tabla_ozone.data_lbl), ...
                        'rowlabels', str2latex(tabla_ozone.events), ...
                        'alignment', 'c');
                    
         % use the aux function changeFigures (see the end of this file) to add the brewer id to the figs' tags
         handles2save=changeFigures(Cal, ["FILTER_DISTRIBUTION_event_","Ozone_diff_filter_rel_comparative_event_"], [], []);
        
         % save the figures using the aux function saveFigs (see the end of this file) to call printfiles_report with
         %+appropriate parameters (in particular, an aux_pattern cell with the file_suffix defined in the input options)
         saveFigs(Cal, handles2save, file_suffix, fig_width, fig_height, fig_linewidth);            
        
        
    end % loop over brewers
    
end % end filters_ozone


%% ******************************************************
%  ******************************************************
%
%  get the ETC filter corrections looking at the langleys
%
%  ******************************************************
%  ******************************************************

if filters_langley_run
    
    % load config, ozone_ds, and ozone_raw data
    % this can be done with a call to the aux function load_bdata or loading a mat file where the data was previously saved
    if filters_langley_load_data
        load(filters_langley_data_file, 'config', 'ozone_ds' , 'ozone_raw');
    else
        [config, ozone_ds, ozone_raw]=load_bdata(Cal, filters_langley_data_file);
    end
    
    %% loop over all the brewers
    for brw_idx=brewers2check %Cal.n_ref
        close all;
            
        % the brewer to analyze is Cal.brw{Cal.n_inst}
        Cal.n_inst=brw_idx;
        fprintf(strcat("\n\n******** running filters_langley for Brewer #",Cal.brw_str{Cal.n_inst},"\n"))
    
        % initialize output table
        tabla_langley.data=[];
        tabla_langley.data_lbl={'F#1 corr', 'F#2 corr', 'F#3 corr', 'F#4 corr', 'F#5 corr'};
        tabla_langley.events=filter_events{brw_idx}.labels;
        
        %% do the langley calculation including the filter data
        % format input data
        [ozone_lgl{Cal.n_inst}, cfg_indv, leg, ozone_lgl_sum{Cal.n_inst}] = langley_data_cell(ozone_raw{Cal.n_inst}, ...
                                                                                              ozone_ds{Cal.n_inst}, ...
                                                                                              config{Cal.n_inst});

        % regression for the half-day langleys, taking into account filter groups 0-1-2, 3, 4, and 5
        % some args are
        %     'O3_hday' -> max std for each half day
        %     'N_hday'  -> min num of summaries for each half day
        %     'F_corr'  -> keep at [0  0  0  0  0  0] to get the real ETCs per filter
        %     'AOD'     -> Aeronet file to filter with the AOD
        %     'Cloud'   -> Cloud screening file to filter clouds
        my_args={ozone_lgl{Cal.n_inst}, 'plots', 0, 'airmass', [1.05 4], 'O3_hday', 2, 'N_hday', 12, ...
                 'F_corr', [0  0  0  0  0  0]};
             
        if ~isempty(filters_langley_aod_file)
            my_args{end+1}='AOD';
            my_args{end+1}=filters_langley_aod_file;
        end
        
        if ~isempty(filters_langley_cloud_file)
            my_args{end+1}='Cloud';
            my_args{end+1}=filters_langley_cloud_file;
        end
        
        lgl_filt{Cal.n_inst} = langley_filter_lvl1(my_args{:});

        % ETC filter corrections
        brw_indv_ = langley_analys_filter(lgl_filt, Cal.n_inst, ...
                                          'res_filt', 1, 'plot_flag', 0);

        % the previous codes return data for two configs, but we will analyze only the operative (should be no. 1)
        cfgs=1;      % operative cfg
        
        % format ouput
        nd0_ = cat(1, brw_indv_(:, [1 2], cfgs), brw_indv_(:, [1 6], cfgs) ); 
        nd0 = sortrows(nd0_, 1);
        nd3_ = cat(1, brw_indv_(:, [1 3], cfgs), brw_indv_(:, [1 7] ,cfgs) ); 
        nd3 = sortrows(nd3_, 1);       
        nd4_ = cat(1, brw_indv_(:, [1 4], cfgs), brw_indv_(:, [1 8], cfgs) );
        nd4 = sortrows(nd4_, 1);

        %% table & histogram with data per event
        num_events=length(filter_events{brw_idx}.dates);
        for event_idx=1:num_events
            fprintf(strcat("**** working on event #",num2str(event_idx),"\n"))
            
            % slice by events
            nd0_slice=sliceVarByEvent(nd0, filter_events{brw_idx}.dates, event_idx, num_events);
            nd3_slice=sliceVarByEvent(nd3, filter_events{brw_idx}.dates, event_idx, num_events);
            nd4_slice=sliceVarByEvent(nd4, filter_events{brw_idx}.dates, event_idx, num_events);
            
            % median values of the half-day ETCs for each filter
            nd0_etc=nanmedian(nd0_slice(:,2));
            nd3_etc=nanmedian(nd3_slice(:,2));
            nd4_etc=nanmedian(nd4_slice(:,2));
        
            % filter corrections are calculated wrt the ETC of nd#0
            nd3_etc_corr=round(nd3_etc-nd0_etc,0); % original code from Langley185.m
            nd4_etc_corr=round(nd4_etc-nd0_etc,0);
            
            %nd3_etc_corr=round( nanmedian( nd3_slice(:,2) - nd0_slice(:,2) ) , 0); % calc using all data 
            %nd4_etc_corr=round( nanmedian( nd4_slice(:,2) - nd0_slice(:,2) ) , 0);
        
            % print data table
            if verbose_output
                fprintf(strcat("** Median of the ETC of each filter over all the Langleys of event ", ...
                                num2str(event_idx), "\n"))
                            
                fprintf(strcat("** ETC for ND #0-1-2 = ", num2str(nd0_etc) ,"\n") )
                fprintf(strcat("** ETC for ND #3     = ", num2str(nd3_etc) ) )
                fprintf(strcat(" -> ETC filter correction = ", num2str(nd3_etc_corr), "\n") )
            
                if length(nd4_slice(~isnan(nd4_slice(:,2)),2))>2 % nhist requires more than 2 datapoints
                    fprintf(strcat("** ETC for ND #4     = ", num2str(nd4_etc) ) )
                    fprintf(strcat(" -> ETC filter correction = ", num2str(nd4_etc_corr), "\n") )
                end
            end
            
            % save data of this event to a table, to save it at the end of the events' loop
            tabla_langley.data(event_idx,1)=filter_events{brw_idx}.dates(event_idx);
            tabla_langley.data(event_idx,2)=NaN; % nothing for filter #1
            tabla_langley.data(event_idx,3)=NaN; % nothing for filter #2
            tabla_langley.data(event_idx,4)=nd3_etc_corr; % filter #3
            if length(nd4_slice(~isnan(nd4_slice(:,2)),2))>2
                tabla_langley.data(event_idx,5)=nd4_etc_corr; % filter #4
            else
                tabla_langley.data(event_idx,5)=NaN; % nothing for filter #4
            end
            tabla_langley.data(event_idx,6)=NaN; % nothing for filter #5
            
            % histogram of ETCs for each filter
            figure; 
            set(gcf,'Tag',strcat('ETC_FILTER_HISTOGRAM_event_',num2str(event_idx)));

            % two cases: with and without data for filter 4
            if length(nd4_slice(~isnan(nd4_slice(:,2)),2))>2
                aux={nd0_slice(~isnan(nd0_slice(:,2)),2), nd3_slice(~isnan(nd3_slice(:,2)),2), nd4_slice(~isnan(nd4_slice(:,2)),2)};
            else
                aux={nd0_slice(~isnan(nd0_slice(:,2)),2), nd3_slice(~isnan(nd3_slice(:,2)),2)};
            end

            nhist(aux,'box','smooth','samebins','ylabel','PDF','xlabel','ETC');
            fprintf("\n") % fix nhist always writing a warning without carriage return!
        
            title(sprintf('%s Filter ETCs (%s to %s)', ...
                  Cal.brw_name{Cal.n_inst}, datestr(nd3_slice(1,1),29), datestr(nd3_slice(end,1),29) ) );
            
            grid; 
            box on;

            set(findobj(gca,'Type','Line'),'Marker','None','LineWidth',2)

            % final touches to the fig, depending on the filter 4 data being present or not
            if length(nd4_slice(~isnan(nd4_slice(:,2)),2))>2
                legend({'ND#0,1,2','ND#3', 'ND#4'}, 'Location', 'Best', 'Orientation', 'Vertical');
                vline_bottom([nd0_etc nd3_etc nd4_etc], '--k', ...
                        {num2str(round(nd0_etc,0)), num2str(round(nd3_etc,0)), num2str(round(nd4_etc,0))} )
            else
                legend({'ND#0,1,2','ND#3'}, 'Location', 'Best', 'Orientation', 'Vertical');
                vline_bottom([nd0_etc nd3_etc], '--k', {num2str(round(nd0_etc,0)), num2str(round(nd3_etc,0))} )
            end       
                        
        end % loop of events
            
        %% change the sign of table_langley to have the same as in table_fioavg
        tabla_langley.data(:,2:end)=-tabla_langley.data(:,2:end);
        
        % show again the Filter ETC corr per event but this time as a table
        fprintf("**** ETC Filter corr from the Langley calculation, taking as reference the ETC of filter group #0,1,2\n")
        displaytable(tabla_langley.data(:,2:end), tabla_langley.data_lbl, 10, '.2f', tabla_langley.events);
        
        % save table
        fprintf(strcat("** saving table to ", Cal.dir_tables, "\n"))
        
        % make sure file_suffix is a vector of chars or the following line will break!
        table_file=fullfile(Cal.dir_tables, ['table_ETC_FILTER_CORR_LANGLEY_', Cal.brw_str{Cal.n_inst}, '_', file_suffix, '.tex'] );
        
        matrix2latex_ctable(tabla_langley.data(:,2:end), table_file, ...
                        'columnlabels', str2latex(tabla_langley.data_lbl), ...
                        'rowlabels', str2latex(tabla_langley.events), ...
                        'alignment', 'c');
         
        % save figs
        handles2save=changeFigures(Cal, "ETC_FILTER_HISTOGRAM_event_", [], []);
        
        saveFigs(Cal, handles2save, file_suffix, fig_width, fig_height, fig_linewidth);
               
        %% Time Series of ETC filter corrections      
        % calc median, standard error of the mean, and std of the ETC filter corrections, for every week of every year
        % (the original code from Langley185.m used the mean instead of the median, but the median is used in the previous section,
        % so I made the change to get more consistent results)
        
        % note i'm changing the signs of the differences to get the same as in the FIOAVG and ozone analyses!
        [m_brw, s_brw, n_brw] = grpstats( [nd0(:,1) -nd3(:,2)+nd0(:,2) -nd4(:,2)+nd0(:,2)], ...
                                          {year(nd0(:,1)) weeknum(nd0(:,1))}, {'mean','sem','std'} );                                      

        % uncomment the next two lines to get the ETC (not ETC corr) plot
        %[m_brw, s_brw, n_brw] = grpstats( [nd0(:,1) nd3(:,2) nd4(:,2)], ...
        %                                  {year(nd0(:,1)) weeknum(nd0(:,1))}, {'median','sem','std'} );  
        
        
        % separate the output of grpstats in three vars
        lmu=sortrows(m_brw,1);
        lsem=sortrows(cat(2,m_brw(:,1),s_brw(:,2:end)),1); 
        lstd=sortrows(cat(2,m_brw(:,1),n_brw(:,2:end)),1);
   
        % init plot
        figure; 
        set(gcf,'Tag','ETC_FILTER_SERIES');
        
        grid;
        
        title(sprintf('%s ND Filter Corr. (%s to %s)', ...
              Cal.brw_name{Cal.n_inst}, datestr(lmu(1,1),29), datestr(lmu(end,1),29) ) );
          
        ylabel('ETC Filter Correction');
          
        % data for filter 3
        idx=lsem(:,2)==0 | isnan(lsem(:,2)); 
        lsem3=lsem(~idx,:);  
        lmu3=lmu(~idx,:);
        h_nd3=boundedline(gca,lmu3(:,1),lmu3(:,2),lsem3(:,2),'-b','alpha' ,'transparency', 0.3);
             
        % data for filter 4
        idx=lsem(:,3)==0 | isnan(lsem(:,3)); 
        lsem4=lsem(~idx,:);  
        lmu4=lmu(~idx,:);

        % check if there is data for filter 4 and complete the figure
        if ~isempty(lmu4)
            h_nd4=boundedline(gca,lmu4(:,1),lmu4(:,3),lsem4(:,3),'-r','alpha' ,'transparency', 0.3);
            
            legend([h_nd3, h_nd4], {'ND#3', 'ND#4'}, 'Location', 'Best', 'Orientation', 'Vertical');
        else
            legend(h_nd3, 'ND#3', 'Location', 'Best', 'Orientation', 'Vertical');
        
        end
           
        % add horizontal lines at the events
        num_events=length(filter_events{brw_idx}.dates);
        for event_idx=1:num_events
            
            lmu3_slice=sliceVarByEvent(lmu3, filter_events{brw_idx}.dates, event_idx, num_events);
            lmu4_slice=sliceVarByEvent(lmu4, filter_events{brw_idx}.dates, event_idx, num_events);   
                        
            nd3_etc_corr_series=round(median(lmu3_slice(:,2)));
            hline_pos( nd3_etc_corr_series, '--b', num2str(nd3_etc_corr_series), 0.9*(event_idx)/(num_events+1) );
            
            nd4_etc_corr_series=round(median(lmu4_slice(:,3)));
            hline_pos( nd4_etc_corr_series, '--r', num2str(nd4_etc_corr_series), 0.9*(event_idx)/(num_events+1) );
        end
           
        datetick('x',12,'KeepTicks');
        box on;
        
        % reduce y axis limits if too large
        my_ylim=ylim;
        if my_ylim(2)>30
            ylim([my_ylim(1),30]);
        end
        
        my_ylim=ylim;
        if my_ylim(1)<-30
            ylim([-30,my_ylim(2)]);
        end
        
        % save figs
        handles2save=changeFigures(Cal, "ETC_FILTER_SERIES", fig_xlim, filter_events{brw_idx});
        
        saveFigs(Cal, handles2save, file_suffix, fig_width, fig_height, fig_linewidth);
               
    end % loop over brewers

end % filters_langley


%% ********************************************************************************
%  ********************************************************************************
%
%                                 auxiliary functions
%
%  MATLAB does not find these functions if they are called from the command window
%  or using F9, but it founds them if the whole program or a section is executed...
%
%  ********************************************************************************
%  ********************************************************************************

function handles2save=changeFigures(Cal, figure_tags, fig_xlim, events)
    % modify figures: add the brewer id to the tag, set the limit of the x axis, add events

    % dimensions of handles2save must be 1 row, multiple cols or printfile_report won't work
    % length(figure_tags) will be the correct number of figs only if there are no multiple figs starting with the same tag
    handles2save=zeros(1,length(figure_tags)); 

    fig_idx=1;    
    for my_tag=figure_tags
        % look for all the figs with tags STARTING by my_tag
        my_tag_regexp=strcat('^',my_tag);
        
        figure_handles=findobj('-regexp','Tag',my_tag_regexp);
    
        % there might be multiple handles starting by the same tag, loop over all of them       
        for my_fig_handle_idx=1:length(figure_handles)
                       
            handles2save(fig_idx)=figure(figure_handles(my_fig_handle_idx)); % activate and save this handle
            
            fig_tag=get(gcf,'tag');
            fig_tag=strcat(Cal.brw_str{Cal.n_inst},'_',fig_tag);
            set(gcf,'Tag',fig_tag);
    
            if ~isempty(fig_xlim)
                xlim(fig_xlim)
            end
        
            if ~isempty(events)
                vline_v(events.dates,'k',events.labels)
            end
            
            fig_idx=fig_idx+1;
        end % multiple handles starting with the same tag
            
    end % tags
    
end

% ----------------------------
function saveFigs(Cal, figure_handles, file_suffix, fig_width, fig_height, fig_linewidth)
    % a printfiles_report wrapper, mostly to add an aux_pattern cell
    
    % some notes on the behaviour of printfiles_report:
    %
    % 1) if the output_dir contains something like "185_figures", it will be used as
    %+the start of the filename used to save the file; otherwise, the filename will
    %+start with the string "General"
    %
    % 2) the tags are also used as part of the filename. files won't be overwritten if
    %+multiple handles with the same tag are passed, because "_1", "_2", ... is added
    %+by default to the filename of the second, third, ... figures
    %
    % 3) you can pass an "aux_pattern" optional argument, which is a cell with a suffix
    %+for each fig_handle passed
    %
    % 4) figure_handles must be a 1 row vector with as many cols as needed (1 col and many rows won't work)
    
        
    % define an aux_pattern to avoid a mess of _1, _2, ... suffixes
    fig_suffix_cell=cell(length(figure_handles),1);
    fig_suffix_cell(:)={file_suffix}; % make sure file_suffix is a vector of chars so that printfiles doesn't break!
        
    fprintf(strcat("** saving figs to ",Cal.dir_figs,"\n"))
       
    %disp(figure_handles)
    
    printfiles_report(figure_handles, Cal.dir_figs, ...
        'Width', fig_width, 'Height', fig_height, 'LineWidth', fig_linewidth, ...
        'aux_pattern', fig_suffix_cell);
end

% -----------------------------
function [summary, summary_old]=load_summ(brw_idx, Cal)
    % load the summaries from files such as files such as CODE/iberonesia/RBCC_E/Triad/Langley/summary_Brw1_20162019.txt
    % note these files are generated by CODE/iberonesia/RBCC_E/Triad/Langley/load_data_2016_2019.m 
    %+(or is it CODE/iberonesia/RBCC_E/Triad/Langley/load_data.m?)
    % note also that CODE/iberonesia/RBCC_E/Triad/Langley/read_summaies_2016_2019.m seems to generate yearly txt files in the yearly dirs

    summ_file_prefix=["summary_Brw", "summary_old_Brw"];
    
    summ_data=cell(2,1);
    
    for summ_idx=1:2
        summ_table=[char(summ_file_prefix(summ_idx)), num2str(brw_idx), '_20162019.txt'];
        summ_table=fullfile(Cal.path_root,'..','Triad','Langley',summ_table);
        
        fprintf(strcat("** loading data from summary file ", summ_table,"\n") )
        
        summ_table=readtable(summ_table);
        
        summ_data{summ_idx}=[summ_table.x_Date, summ_table.sza, summ_table.m2, summ_table.temp, summ_table.nd, ...
                             summ_table.O3_1, summ_table.std, summ_table.ms9_corr, summ_table.ms9, summ_table.O3_2, ...
                             summ_table.std, summ_table.O3_1_sl, summ_table.std, summ_table.R6_ref, ...
                             summ_table.R6_calc, summ_table.A1, summ_table.ETC];
    end
    
    summary=summ_data{1};
    summary_old=summ_data{2};    
end

% -----------------------------
function [config, ozone_ds, ozone_raw]=load_bdata(Cal, file2save)
    % load data from B files
    % adapted from section "READ Brewer Summaries" in CODE/iberonesia/RBCC_E/Triad/Langley/load_data.m
    
    for brw_idx=Cal.n_ref
        fprintf(strcat("** loading B files from Brewer #", Cal.brw_str{Cal.n_inst},"\n"))
        %ozone_sum{i}={};
        config{brw_idx}={};
        ozone_ds{brw_idx}={}; 
        ozone_raw{brw_idx}={};
        %ozone_raw0{i}={};
        %sl{i}={};
        %sl_cr{i}={};

        [ozone,log_,missing_]=read_bdata(brw_idx,Cal);
    
        % depuramos datos (ver incidencias en config. matrix)
        %ozone=dep_data(Cal.incidences_text{i},ozone);

        %ozone_sum{i}=ozone.ozone_sum;
        config{brw_idx}=ozone.config;
        ozone_ds{brw_idx}=ozone.ozone_ds;
        ozone_raw{brw_idx}=ozone.raw;
        %ozone_raw0{i}=ozone.raw0;
        %sl{i}=ozone.sl;       % first calibration / bfiles
        %sl_cr{i}=ozone.sl_cr; % recalc. with 2? configuration  
        
        fprintf(strcat("** saving variables config, ozone_ds, and ozone_raw to file ", file2save))
        save(file2save, 'config', 'ozone_ds', 'ozone_raw');
    end    
end

% ---------------------------------
function hhh=vline_bottom(x,in1,in2)
    % vline_v modified to print labels close to the bottom 
    % vline_v by Brandon Kuczenski for Kensington Labs, brandon_kuczenski@kensingtonlabs.com

    if length(x)>1  % vector input
        for I=1:length(x)
            switch nargin
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
            end
            h(I)=vline_bottom(x(I),linetype,label);
        end
    else
        switch nargin
        case 1
            linetype='r:';
            label=num2str(x);
        case 2
            linetype=in1;
             label=num2str(x);
        case 3
            linetype=in1;
            label=in2;
        end  
    
        g=ishold(gca);
        hold on

        y=get(gca,'ylim');
        h=plot([x x],y,linetype);
        if length(label)
            xx=get(gca,'xlim');
            xrange=xx(2)-xx(1);
            xunit=(x-xx(1))/xrange;
            if xunit<0.8
                text(x,y(1)+0.2*(y(2)-y(1)),label,'color',get(h,'color'),'Rotation',-90,'BackgroundColor','w','FontWeight','Bold');
            else
                text(x,y(1)+0.2*(y(2)-y(1)),label,'color',get(h,'color'),'Rotation',-90,'BackgroundColor','w','FontWeight','Bold');
            end
        end     

        if g==0
        hold off
        end
        set(h,'tag','vline','handlevisibility','off')
    end % else

    if nargout
        hhh=h;
    end

end

% ----------------------------------------
function hhh=hline_pos(y,in1,in2,pos)
    % hline modified to put the label on top of the line and at x position specified by pos (0.1 ~ left, 0.5 ~ center, 0.9 ~ right) 
    % hline by Brandon Kuczenski for Kensington Labs, brandon_kuczenski@kensingtonlabs.com

    if length(y)>1  % vector input
        for I=1:length(y)
            switch nargin-1
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
            end
            h(I)=hline_pos(y(I),linetype,label,pos);
        end
    else
        switch nargin-1
        case 1
            linetype='r:';
            label='';
        case 2
            linetype=in1;
            label='';
        case 3
            linetype=in1;
            label=in2;
        end 
    
        g=ishold(gca);
        hold on
        
        x=get(gca,'xlim');
        h=plot(x,[y y],linetype);
        if ~isempty(label)
            yy=get(gca,'ylim');
            yrange=yy(2)-yy(1);
            yunit=(y-yy(1))/yrange;
            if yunit<0.2            
                text(x(1)+pos*(x(2)-x(1)),y,label,'color',get(h,'color'),'BackgroundColor','w')
            else
                text(x(1)+pos*(x(2)-x(1)),y,label,'color',get(h,'color'),'BackgroundColor','w')
            end
        end

        if g==0
            hold off
        end
        set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
    end % else

    if nargout
        hhh=h;
    end
    
end


% -------------------------------
function slice=sliceVarByEvent(var, event_dates, event_num, total_num_events)
    % slice (subset) a variable according to a set of events

    if event_num ~= total_num_events
        slice_ini=var(:,1)>=event_dates(event_num);
        slice_end=var(:,1)<=event_dates(event_num+1);      
        slice_idx=slice_ini & slice_end;
    else
        slice_idx=var(:,1)>=event_dates(event_num);
    end
    
    slice=var(slice_idx,:);
    
end


% ---------------------------------
function [results_rel, results_abs]=ozone_filter_analysis_mi_ms9(summary,Cal,n_inst,varargin)
    % this is mostly the same as /CODE/iberonesia/matlab/brewer/ozone/ozone_filter_analysis_mi
    % but returns differences of MS9 ratios, instead of ozone, and it's easier to select the
    % time interval to keep filter changes
    
    %max_time_between_two_measurements=10; % this is in minutes
    max_airmass_between_two_measurements=0.03;
    
    
    %  Analiza el cambio en el ozono durante el cambio de filtros
    %  input variable summary
    %  output: df
    %   medidas consecutivas con cambio de filtro
    %   df=[fecha, filtro_1, filtro_2, summary(filtro_1)-summary(filtro_2)]
    %  
    % Juanjo 14/09/2011: Se a�ade condicional del tipo if isempty() return 
    %                    para salir en caso de no data

    if nargin==2
        summary=summary{Cal.n_inst};
    else
        summary=summary{n_inst};
        Cal.n_inst=n_inst;
    end

    if all(isnan(summary))
        disp('No data');
        return
    end

    aux=sortrows(summary,1); freq_filter=tabulate(aux(:,5));

    cf=diff(aux(:,5));cf=[1;cf];  
    cf_idx=find(cf);% buscamos el cambio de filtro (no nulos) 
    aux=aux(cf_idx,:); cf=cf(cf_idx);
    cf_idx=find(abs(cf)==64);% s�lo filtros consecutivos
    %                 indx     primer F#     segundo F#
    
    %J This is to remove pairs of measurements taking into account the time difference between them
    %chg_filter=[cf_idx,aux(cf_idx-1,5),aux(cf_idx,5),(aux(cf_idx,1)-aux(cf_idx-1,1))*24*60];
    
    %J This is to remove pairs of measurements taking into account the airmass difference between them
    chg_filter=[cf_idx,aux(cf_idx-1,5),aux(cf_idx,5),(aux(cf_idx,3)-aux(cf_idx-1,3))];
    
    %J Keep only data from measurements taken within max_time_between_two_measurements minutes or less
    %chg_filter=chg_filter(chg_filter(:,end)<=max_time_between_two_measurements,:);
    
    %J Keep only data from measurements taken within max_airmass_between_two_measurements airmass or less
    chg_filter=chg_filter(abs(chg_filter(:,end))<=max_airmass_between_two_measurements,:);   
    
    chg_rel=repmat({NaN*ones(fix(size(summary,1)/2),2)},1,8);
    chg_abs=repmat({NaN*ones(fix(size(summary,1)/2),2)},1,8);
 
    % DEFINICIONES DE CAMBIO. No considero cambios de m�s de un filtro 
    % #0 -> #1 (#1 -> #0), #1 -> #2 (#2 -> #1), #2 -> #3 (#3 -> #2), #3 -> #4 (#4 -> #3) 

    id=1;
    for filt=1:4
        % ADELANTE
        idx=chg_filter(chg_filter(:,2)==(filt-1)*64 & chg_filter(:,3)==filt*64,1); idx=sort([idx-1;idx]);
        g=aux(idx,:); primero=g(1:2:end,:);  segundo=g(2:2:end,:); 
        % diff. relativa
        dif_rel=(segundo-primero)*100./primero; dif_rel(:,1)=nanmean([primero(:,1),segundo(:,1)],2);
        chg_rel{id}(1:size(dif_rel,1),:)=dif_rel(:,[1 9]); %J 6->9
        % diff. absoluta
        dif_=segundo-primero; dif_(:,1)=nanmean([primero(:,1),segundo(:,1)],2);
        chg_abs{id}(1:size(dif_,1),:)=dif_(:,[1 9]); %J 6->9

        % ATRAS 
        idx=chg_filter(chg_filter(:,2)==filt*64 & chg_filter(:,3)==(filt-1)*64,1); idx=sort([idx-1;idx]);
        g=aux(idx,:); primero=g(1:2:end,:);  segundo=g(2:2:end,:); 
        % diff. relativa
        dif_rel=(primero-segundo)*100./segundo; dif_rel(:,1)=nanmean([primero(:,1),segundo(:,1)],2);
        chg_rel{id+1}(1:size(dif_rel,1),:)=dif_rel(:,[1 9]); %J 6->9 
        % diff. absoluta    
        dif_=primero-segundo; dif_(:,1)=nanmean([primero(:,1),segundo(:,1)],2); 
        chg_abs{id+1}(1:size(dif_,1),:)=dif_(:,[1 9]); %J 6->9
    
        id=id+2;
    end
    tableform({'F#0 <-> F#1','F#1 <-> F#2','F#2 <-> F#3','F#3 <-> F#4'},...
              [length(find(~isnan(chg_rel{1}(:,1))))+length(find(~isnan(chg_rel{2}(:,1)))),...
               length(find(~isnan(chg_rel{3}(:,1))))+length(find(~isnan(chg_rel{4}(:,1)))),...
               length(find(~isnan(chg_rel{5}(:,1))))+length(find(~isnan(chg_rel{6}(:,1)))),...
               length(find(~isnan(chg_rel{7}(:,1))))+length(find(~isnan(chg_rel{8}(:,1))))]);

    results_rel=cell2mat(chg_rel); results_abs=cell2mat(chg_abs); 
    if isempty(varargin)
       fplot=1; qplot=1;
    elseif size(varargin)==1
       fplot=varargin{1};   qplot=0;
    else
       fplot=varargin{1};   qplot=varargin{2};   
    end

       % queso con porcentajes relativos al total
       if qplot
          figure; set(gcf,'Tag','FILTER_DISTRIBUTION');
          label_1=mmcellstr(sprintf('F#%d= |',freq_filter(:,1)./64));
          label_=mmcellstr(sprintf('%.1f%% |',freq_filter(:,3)));
          explode = repmat(1,1,size(freq_filter,1))'; pie3(freq_filter(:,2),explode,strcat(label_1,label_));
          set(findobj(gcf,'Type','text'),'Backgroundcolor','w');
          title(sprintf('%s%s','Filters:  ',Cal.brw_name{Cal.n_inst}),'FontSize',12,'FontWeight','Bold');
       end
   
    if fplot
       f=figure;  set(f,'Tag','Ozone_diff_filter_rel');
       rectangle('Position',[.5,-0.5,8,1],'FaceColor',[.95 .95 .95]); hold on; 
       boxplot(results_rel(:,2:2:end),'notch','on',...% ploteamos campo 6
                        'labels',{'0-64','64-0','64-128','128-64','128-192','192-128','192-256','256-192'}); 
       set(gca,'YLim',[-2 2]); ylabel('Relative Difference (%)');  hline(0,'-.k');  grid;
       title(sprintf('Ozone difference by filter chg. %s\r\n Referenced always to lower filter for each group',Cal.brw_name{Cal.n_inst}));

       f=figure;  set(f,'Tag','Ozone_diff_filter');   hold on; 
       boxplot(results_abs(:,2:2:end),'notch','on',...% ploteamos campo 6
                         'labels',{'0-64','64-0','64-128','128-64','128-192','192-128','192-256','256-192'}); 
       ylabel('Absolute Difference');  hline(0,'-.k');  grid;
       title(sprintf('Ozone difference by filter chg. %s',Cal.brw_name{Cal.n_inst})); 
    end
    
end