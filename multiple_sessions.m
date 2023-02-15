%% clear workspace
clear all; close all; clc

%% list the desired sessions
session_names = {'2021_6_25_2000_BehaviorToneWater.mat',
    '2021_7_1_2130_BehaviorToneWater.mat',
    '2021_7_10_2300_BehaviorToneWater.mat'
    '2021_7_18_2300_BehaviorToneWater.mat',
    '2021_7_25_2300_BehaviorToneWater.mat',
    '2021_8_7_2300_BehaviorToneWater.mat'};

%% load the saved data
n_sessions = numel(session_names);
session_data = cell(n_sessions,1);
lead_lag_vec = -1000:20:1000;

folder = 'C:\Users\selim\Dropbox\dannce-analysis-code\cells_examples\sfn_figures\';

%%
counter_1 = 0;
for a=1:n_sessions
   session_data{a} = load(session_names{a});
   mkdir(folder, session_names{a});
   total_cell_number = numel(session_data{a}.cell_type);
   counter_1 =  total_cell_number + counter_1;
end

total_data = zeros(counter_1,6);

counter = 0;
for s=1:n_sessions
    saved_folder = [folder session_names{s} '\']
    num_cells = numel(session_data{s}.cell_type);
    
    %1)lead-lag business for individual cells
    for k = 1:num_cells
        smooth_fit = smoothdata(session_data{s}.avg_test_fit(:,k),'gaussian',10);
        figure(1)
        plot(lead_lag_vec,smooth_fit,'k')
        hold on
        plot([0 0],[min(smooth_fit) max(smooth_fit)],'-r')
        title(k)
        xlabel('Time shift')
        ylabel('Pseudo r^2 values') 
        hold off
        %pause 
        box off
        baseFileName = sprintf('Cell%d', k); %if you want to save figures
        %saveas(gcf, [saved_folder baseFileName],'pdf');
    end
    
    %3)histogram of cells not fit vs. fit
    [max_r2_value, max_r2_index] = max(session_data{s}.avg_test_fit, [], 1);
    [cell_type_s,peak_width_s,mean_firing_rate_s,peak_to_valley_s] = classify_cells(session_data{s}.ephys_struct_pr);
 
    total_data(counter+1:counter+num_cells,1) = max_r2_value;
    total_data(counter+1:counter+num_cells,2) = max_r2_index;
    total_data(counter+1:counter+num_cells,3) = cell_type_s;
    total_data(counter+1:counter+num_cells,4) = peak_width_s;
    total_data(counter+1:counter+num_cells,5) = mean_firing_rate_s;
    total_data(counter+1:counter+num_cells,6) = peak_to_valley_s;
    
    min_max_r2 = min(max_r2_value);
    max_max_r2 = max(max_r2_value);
    [n,x] = hist(max_r2_value, linspace(min_max_r2,max_max_r2,10));
    bar(x,n,'k')
    box off
    set(gca,'fontsize',20)
    axis([min_max_r2-0.1 max_max_r2+0.1 0 max(n)+2])
    hold on
    plot([0 0],[0 max(n)],'-r', 'LineWidth', 3)
    xlabel('Max Pseudo r^2 values')
    ylabel('# of cells')
    title('Total cell fit')
    saveas(gcf, [saved_folder 'Model_success_hist'],'pdf');
    hold off
    
    %4) histogram of leading vs. lagging cells
    %restrict it to 'good cells', i.e. cells that have positive pseudo r^2
    good_cells_vec = max_r2_value > 0; %'good cell'
    num_cells_good = sum(good_cells_vec);
    ind_of_cells_good = find(good_cells_vec);
    avg_test_fit_good = session_data{s}.avg_test_fit(:,good_cells_vec);
    [max_r2_good, max_r2_ind_good] = max(avg_test_fit_good,[],1);
    
    [n_num,x_shift] = hist(lead_lag_vec(max_r2_ind_good),linspace(min(lead_lag_vec),max(lead_lag_vec),20));
    bar(x_shift,n_num,'k')
    box off
    set(gca,'fontsize',20)
    axis([-1050 inf 0 max(n_num)+2])
    hold on
    plot([0 0],[0 max(n_num)],'r','linewidth',2)
    xlabel('Time shift')
    ylabel('# of cells')
    title('Leading vs. lagging cells')
    saveas(gcf, [saved_folder 'lead_vs_lag_cells'],'pdf');
    hold off
    
    %'heatmap' thingy
    start_end_dannce_attempt_s = session_data{s}.start_end_dannce_attempt;
    attempt_length = start_end_dannce_attempt_s(:,2)-start_end_dannce_attempt_s(:,1);
    max_length = max(attempt_length);
    [sorted_attempts, ind_att] = sort(attempt_length);
    n_attempts = size(start_end_dannce_attempt_s,1);
    avg_attempt_len = round(mean(attempt_length));
    endpoint = max_length + 50; 
    image_data = zeros(num_cells_good, avg_attempt_len-1); 

    for num = 1:num_cells_good
        st_num = session_data{s}.ephys_struct_pr(ind_of_cells_good(num)).st;
        trial_fr = nan(n_attempts,avg_attempt_len-1);
        for a = 1:n_attempts
            new_ts_trial = linspace(session_data{s}.t_us(start_end_dannce_attempt_s(a,1)+(lead_lag_vec(max_r2_ind_good(num))/10)),session_data{s}.t_us(start_end_dannce_attempt_s(a,2)),avg_attempt_len+2);
            [n_h,x_h] = hist(st_num, new_ts_trial);
            s_per_bin = median(diff(x_h))/1000;
            rebinned_at = n_h(2:avg_attempt_len); % re-binned spiketrain
            trial_fr(a,:) = rebinned_at./s_per_bin;
        end
        image_data(num,:) = smoothdata(mean(trial_fr, 1),'gaussian',10);
        close all
    end

    % plot heatmap thingy
    data_z = zscore(image_data')';
    [val,ind] = max(data_z,[],2);
    [val2,ind2]  = sort(ind,'ascend');
    imagesc(data_z(ind2,:))
    xlabel('Warped time');
    ylabel('The cell #')
    colorbar
    saveas(gcf, [saved_folder 'firing_rate_time'],'pdf');
    hold off
    
    % max shift value vs. time durign attempt
    [tuning_curve,occupancy,st_dev,xvec] = compute_1d_tuning_curve(ind,lead_lag_vec(max_r2_ind_good),6,0,135);
    errorbar(xvec,tuning_curve,st_dev./sqrt(occupancy),'k','linewidth',1.5)
    hold on
    plot([0 135],[0 0],'--b')
    xlabel('Warped time');
    ylabel('The time shift')
    hold off
    saveas(gcf, [saved_folder 'max_shift_time'],'pdf');
    
    
    for num_r = 1:num_cells_good
        st_num = session_data{s}.ephys_struct_pr(ind_of_cells_good(num_r)).st;
        [n_r,x_r] = hist(st_num, session_data{s}.t_us);
        for a = 1:n_attempts
            endpoint_plot = endpoint - sorted_attempts(a);
            spike_time = find(n_r(start_end_dannce_attempt_s(ind_att(a),1)-50+(lead_lag_vec(max_r2_ind_good(num_r))/10):start_end_dannce_attempt_s(ind_att(a),2)+endpoint_plot+(lead_lag_vec(max_r2_ind_good(num))/10))>0);
            end_attempt = sorted_attempts(a)+51; 
            plot(spike_time, a*ones(size(spike_time)),'.k')
            hold on
            plot([end_attempt end_attempt], [a-1 a], '-b')
        end
        plot([50 50], [1 n_attempts] ,'-r')
        xlabel('Binned time')
        ylabel('Attempt #')
        title('Shifted Spike trains of the Cell')
        hold off
        baseFileName2 = sprintf('Spiking_Cell%d', num_r);
        %saveas(gcf, [saved_folder baseFileName2],'pdf');
        close all
    end

    
    %find cell types
    ind_ephys_pr_cells = 1:1:num_cells;
    spn_cell_ind = ind_ephys_pr_cells(session_data{s}.cell_type == 0);
    fsi_cell_ind = ind_ephys_pr_cells(session_data{s}.cell_type == 1);
    good_fsi_cell_ind = intersect(ind_of_cells_good, fsi_cell_ind);
    good_spn_cell_ind = intersect(ind_of_cells_good, spn_cell_ind);
    lead_lag_fsi = lead_lag_vec(max_r2_index(good_fsi_cell_ind));
    lead_lag_spn = lead_lag_vec(max_r2_index(good_spn_cell_ind));
    
    %spn cells - how well fit they are by the model
    max_r2_value_spn = max_r2_value(spn_cell_ind);
    min_max_r2_spn = min(max_r2_value_spn);
    max_max_r2_spn = max(max_r2_value_spn);
    [n_spn,x_spn] = hist(max_r2_value_spn, linspace(min_max_r2_spn, max_max_r2_spn,20));
    bar(x_spn,n_spn,'k')
    box off
    set(gca,'fontsize',20)
    axis([min_max_r2_spn-0.1 max_max_r2_spn+0.1 0 max(n_spn)+1])
    hold on
    plot([0 0],[0 max(n_spn)],'-r', 'LineWidth', 3)
    xlabel('Max Pseudo r^2 values')
    ylabel('# of cells')
    title('Cell fit for spns')
    saveas(gcf, [saved_folder 'Spn_Model_success_hist'],'pdf');
    hold off
    
    %spn cells hist
    [n_num_spn,x_shift_spn] = hist(lead_lag_spn,linspace(min(lead_lag_spn),max(lead_lag_spn),20));
    bar(x_shift_spn,n_num_spn,'k')
    box off
    set(gca,'fontsize',20)
    axis([-1050 1000 0 max(n_num_spn)+2])
    hold on
    plot([0 0],[0 max(n_num_spn)],'r','linewidth',2)
    xlabel('Time shift')
    ylabel('# of cells')
    title('Leading vs. lagging cells for spn')
    saveas(gcf, [saved_folder 'spn_cells'],'pdf');
    hold off
    
    %fsi cells how well fit they are by the model
    max_r2_value_fsi = max_r2_value(fsi_cell_ind);
    min_max_r2_fsi = min(max_r2_value_fsi);
    max_max_r2_fsi = max(max_r2_value_fsi);
    [n_fsi,x_fsi] = hist(max_r2_value_fsi, linspace(min_max_r2_fsi, max_max_r2_fsi,10));
    bar(x_fsi,n_fsi,'k')
    box off
    set(gca,'fontsize',20)
    axis([min_max_r2_fsi-0.1 max_max_r2_fsi+0.1 0 max(n_fsi)+1])
    hold on
    plot([0 0],[0 max(n_fsi)],'-r', 'LineWidth', 3)
    xlabel('Max Pseudo r^2 values')
    ylabel('# of cells')
    title('Cell fit for fsi')
    saveas(gcf, [saved_folder 'Fsi_Model_success_hist'],'pdf')
    hold off
    
    %fsi cells hist
    [n_num_fsi,x_shift_fsi] = hist(lead_lag_fsi,linspace(min(lead_lag_fsi),max(lead_lag_fsi),10));
    bar(x_shift_fsi,n_num_fsi,'k')
    box off
    set(gca,'fontsize',20)
    axis([-1050 1000 0 max(n_num_fsi)+2])
    hold on
    plot([0 0],[0 max(n_num_fsi)],'r','linewidth',2)
    xlabel('Time shift')
    ylabel('# of cells')
    title('Leading vs. lagging cells for fsi')
    saveas(gcf, [saved_folder 'fsi_cells'],'pdf');
    hold off
    
    %fsi when they spike during behavior
    
    
    %spn when they spike during behavior
    
    

    
    counter =  num_cells + counter;
end
%% total data
min_max_r2_total = min(total_data(:,1));
max_max_r2_total = max(total_data(:,1));
[n_total_fit,x_total_fit] = hist(total_data(:,1), linspace(min_max_r2_total,max_max_r2_total,10));
bar(x_total_fit,n_total_fit,'k')
box off
set(gca,'fontsize',20)
axis([min_max_r2_total-0.1 max_max_r2_total+0.1 0 max(n_total_fit)+2])
hold on
plot([0 0],[0 max(n_total_fit)],'-r', 'LineWidth', 3)
xlabel('Max Pseudo r^2 values')
ylabel('# of cells')
title('Total cell fit')
saveas(gcf, [folder 'Total model fit'],'pdf');
hold off

good_cells_vec_total = total_data(:,1) > 0; %'good cell'
num_cells_good_total = sum(good_cells_vec_total);
ind_of_cells_good_total = find(good_cells_vec_total);
max_r2_good_total = total_data(good_cells_vec_total,1);
max_r2_ind_good_total = total_data(good_cells_vec_total,2);

[n_num_total,x_shift_total] = hist(lead_lag_vec(max_r2_ind_good_total),linspace(min(lead_lag_vec),max(lead_lag_vec),10));
bar(x_shift_total,n_num_total,'k')
box off
set(gca,'fontsize',10)
axis([-1050 inf 0 max(n_num_total)+2])
hold on
plot([0 0],[0 max(n_num_total)],'r','linewidth',2)
xlabel('Time shift')
ylabel('# of cells')
title('Leading vs. lagging cells')
saveas(gcf, [folder 'Total_lead_vs_lag'],'pdf');
hold off

%find cell types
max_r2_value_total = total_data(:,1);
max_r2_shift_total = total_data(:,2);
ind_ephys_pr_cells_total = 1:1:numel(max_r2_value_total);
spn_cell_ind_total = ind_ephys_pr_cells_total(total_data(:,3) == 0);
fsi_cell_ind_total = ind_ephys_pr_cells_total(total_data(:,3) == 1);
good_fsi_cell_ind_total = intersect(ind_of_cells_good_total, fsi_cell_ind_total);
good_spn_cell_ind_total = intersect(ind_of_cells_good_total, spn_cell_ind_total);
lead_lag_fsi_total = lead_lag_vec(max_r2_shift_total(good_fsi_cell_ind_total));
lead_lag_spn_total = lead_lag_vec(max_r2_shift_total(good_spn_cell_ind_total));

%spn cells total - how well fit they are by the model
max_r2_value_spn_total = max_r2_value_total(spn_cell_ind_total);
min_max_r2_spn_total = min(max_r2_value_spn_total);
max_max_r2_spn_total = max(max_r2_value_spn_total);
[n_spn_total,x_spn_total] = hist(max_r2_value_spn_total, linspace(min_max_r2_spn_total, max_max_r2_spn_total,10));
bar(x_spn_total,n_spn_total,'k')
box off
set(gca,'fontsize',20)
axis([min_max_r2_spn_total-0.1 max_max_r2_spn_total+0.1 0 max(n_spn_total)+1])
hold on
plot([0 0],[0 max(n_spn_total)],'-r', 'LineWidth', 3)
xlabel('Max Pseudo r^2 values')
ylabel('# of cells')
title('Cell fit for spns')
saveas(gcf, [folder 'Spn total fit'],'pdf');
hold off

%fsi cells total how well fit they are by the model
max_r2_value_fsi_total = max_r2_value_total(fsi_cell_ind_total);
min_max_r2_fsi_total = min(max_r2_value_fsi_total);
max_max_r2_fsi_total = max(max_r2_value_fsi_total);
[n_fsi_total,x_fsi_total] = hist(max_r2_value_fsi_total, linspace(min_max_r2_fsi_total, max_max_r2_fsi_total,6));
bar(x_fsi_total,n_fsi_total,'k')
box off
set(gca,'fontsize',20)
axis([min_max_r2_fsi_total-0.1 max_max_r2_fsi_total+0.1 0 max(n_fsi_total)+1])
hold on
plot([0 0],[0 max(n_fsi_total)],'-r', 'LineWidth', 3)
xlabel('Max Pseudo r^2 values')
ylabel('# of cells')
title('Cell fit for fsi')
saveas(gcf, [folder 'Fsi total fit'],'pdf');
hold off

%fsi cells hist
[n_num_fsi_total,x_shift_fsi_total] = hist(lead_lag_fsi_total,linspace(min(lead_lag_fsi_total),max(lead_lag_fsi_total),6));
bar(x_shift_fsi_total,n_num_fsi_total,'k')
box off
set(gca,'fontsize',20)
axis([-1050 1000 0 max(n_num_fsi_total)+2])
hold on
plot([0 0],[0 max(n_num_fsi_total)],'r','linewidth',2)
xlabel('Time shift')
ylabel('# of cells')
title('Leading vs. lagging cells for fsi')
saveas(gcf, [folder 'Fsi_lead_lag'],'pdf');
hold off

%fsi when they spike during behavior


%spn when they spike during behavior


%spn cells hist
[n_num_spn_total,x_shift_spn_total] = hist(lead_lag_spn_total,linspace(min(lead_lag_spn_total),max(lead_lag_spn_total),6));
bar(x_shift_spn_total,n_num_spn_total,'k')
box off
set(gca,'fontsize',20)
axis([-1050 1000 0 max(n_num_spn_total)+2])
hold on
plot([0 0],[0 max(n_num_spn_total)],'r','linewidth',2)
xlabel('Time shift')
ylabel('# of cells')
title('Leading vs. lagging cells for spn')
saveas(gcf, [folder 'Spn_lead_lag'],'pdf');
hold off

%% computing the averages
%model fit average for different cell types
mean_fit_fsi_total = nanmean(max_r2_value_fsi_total);
mean_fit_spn_total = nanmean(max_r2_value_spn_total);

%lead-lag average for different cell types
mean_lead_lag_fsi_total = mean(lead_lag_fsi_total);
mean_lead_lag_fsi_total = mean(lead_lag_spn_total);

%%
[n_num_fsi_total,x_shift_fsi_total] = hist(lead_lag_fsi_total,linspace(-1000,1000,15));
[n_num_spn_total,x_shift_spn_total] = hist(lead_lag_spn_total,linspace(-1000,1000,20));
smoothfit_spn = smoothdata(n_num_spn_total, 'gaussian',5);
smoothfit_fsi = smoothdata(n_num_fsi_total, 'gaussian',5);
total_spn = sum(n_num_spn_total);
total_fsi = sum(n_num_fsi_total);
proportion_spn = smoothfit_spn/total_spn;
proportion_fsi = smoothfit_fsi/total_fsi;
plot(x_shift_spn_total,proportion_spn, 'r')
hold on
plot(nanmean(lead_lag_fsi_total)*ones(1,2),[0 max(proportion_fsi)],'--b')
plot(nanmean(lead_lag_spn_total)*ones(1,2),[0 max(proportion_spn)],'--r')
plot(x_shift_fsi_total,proportion_fsi, 'b')
hold off

%%

[n_num_fsi_total,x_shift_fsi_total] = hist(max_r2_value_fsi_total,linspace(-0.4,0.8,20));
[n_num_spn_total,x_shift_spn_total] = hist(max_r2_value_spn_total,linspace(-0.4,0.8,15));
smoothfit_spn = smoothdata(n_num_spn_total, 'gaussian',5);
smoothfit_fsi = smoothdata(n_num_fsi_total, 'gaussian',5);
total_spn = sum(n_num_spn_total);
total_fsi = sum(n_num_fsi_total);
proportion_spn = smoothfit_spn/total_spn;
proportion_fsi = smoothfit_fsi/total_fsi;
plot(x_shift_spn_total,proportion_spn, 'r')
hold on
plot(nanmean(max_r2_value_fsi_total)*ones(1,2),[0 max(proportion_fsi)],'--b')
plot(nanmean(max_r2_value_spn_total)*ones(1,2),[0 max(proportion_spn)],'--r')
plot(x_shift_fsi_total,proportion_fsi, 'b')
hold off
box off
set(gca,'fontsize',20)

%% 3d scatter plot
scatter3(total_data(:,4),total_data(:,5),total_data(:,6),[],total_data(:,3));
saveas(gcf, [folder 'cell_type_scatter'],'pdf');

%% plot the waveforms
x = (0:63)*50; % conversion to us
x_interp = linspace(x(1),x(end),500);
for c=1:num_cells
    M = min(session_data{s}.ephys_struct_pr(c).waveforms, [], 2);

    %get the waveform and interpolate
    [min_val, i_min_row] = min(M);
    waveform_c = ephys_struct_pr(c).waveforms(i_min_row,:);
    waveform_c_interp = interp1(x,waveform_c,x_interp);
    set(gcf,'renderer','painters')
    plot(waveform_c_interp)
    pause

    %{
    get the min index of the waveform
    [min_val, i_min_column] = min(waveform_c_interp);
    [min_val, i_min_row] = min(M);
    waveform_c = ephys_struct(c).waveforms(i_min_row,:);
    waveform_c_interp = interp1(x,waveform_c,x_interp);
    %}
end