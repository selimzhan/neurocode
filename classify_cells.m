function [cell_type_pre,peak_width,mean_firing_rate,peak_to_valley] = classify_cells(ephys_struct)

num_cells = numel(ephys_struct);
cell_type_pre = zeros(1,num_cells); %fsi - 1, spn - 0, unclear - 2, tan - 3

peak_width = nan(num_cells,1);
mean_firing_rate = nan(num_cells,1);
peak_to_valley = nan(num_cells,1);
x = (0:63)*50; % conversion to us
x_interp = linspace(x(1),x(end),500);
x_bin_us = x_interp(2) - x_interp(1);
for c = 1:num_cells
    
    %compute minimum point
    M = min(ephys_struct(c).waveforms, [], 2);
    
    %get the waveform and interpolate
    [min_val, i_min_row] = min(M);
    waveform_c = ephys_struct(c).waveforms(i_min_row,:);
    waveform_c_interp = interp1(x,waveform_c,x_interp);
    
    % get the min index of the waveform
    [min_val, i_min_column] = min(waveform_c_interp);
    
    %compute peak width
    [max_val, peak_to_valley_bin] = max(waveform_c_interp(i_min_column:end));
    peak_to_valley(c) = peak_to_valley_bin*x_bin_us;
    
    %find half min on both sides of the min point
    % baseline = first 20/64 bins of the waveform
    range_c = mean(waveform_c(1:19)) - min_val; range_half = range_c./2;
    half_min = min_val + range_half;
    [cl_half_min_left, cl_ind_min_left] = min(abs(half_min - waveform_c_interp(1:i_min_column)));
    [cl_half_v, cl_index_right] = min(abs(half_min - waveform_c_interp(i_min_column:end)));
    cl_ind_min_right = cl_index_right+i_min_column;
    peak_width(c) = (cl_ind_min_right - cl_ind_min_left)*x_bin_us;
    
    %compute mean_firing rate
    mean_firing_rate(c) = (numel(ephys_struct(c).st)/(ephys_struct(c).st(end)-ephys_struct(c).st(1)))*1000;
    
    
    if peak_width(c) <= 150 && peak_to_valley(c) <= 500 && mean_firing_rate(c) > 0.1
        cell_type_pre(c) = 1;
    elseif peak_width(c) > 150 && peak_to_valley(c) > 500 && mean_firing_rate(c) < 10 
        cell_type_pre(c) = 0;
    else
        cell_type_pre(c) = 2;
    end
    
    if abs(max(waveform_c)) > abs(min(waveform_c))
        cell_type_pre(c) = 2;
    end
    
    %{
    for w=1:4
        figure(1)
        plot(ephys_struct(c).waveforms(w,:));
        hold on
    end
    title(['width = ',num2str(peak_width(c)),' peak2valleytime = ',num2str(peak_to_valley(c)),' cell type = ',num2str(cell_type_pre(c))])
    hold off
    keyboard
    %}
    
    
end

%cell_type = logical(cell_type_pre);

return