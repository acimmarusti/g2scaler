clear; close all
hold off

%range of the data to be compared (in microseconds)
tfin = 4.7;
delay = 0.2;

%set to 1 if you want the program to find the pulse
%If set to zero, then specify valid tini (in microseconds)
findpulse = 1;
tini = 3.5;

%how many different scales to check
ny = 20000;

%Y-scaling range
scy = linspace(1, 0.3, ny);

%Y offset options:
% 0 - NO offset
% 1 - Mean offset - 1.0
% 2 - Curve offset
option = 2;

%Size of bins in data (default: 100)
binsize = 100;

%Minumum time resolution (in microseconds)
min_time_res = 164.61e-6;

% Fraction of distance between peaks 
% 0 < frac <= 1
% 1 -> full distance between peaks
min_peak_frac = 0.7;

%total fitting parameters:
%yscaling, xshift, upper envelope(3) and lower envelope(3)
fit_param = 8;

min_time_resolve = binsize * min_time_res;

hold on

if (ispc)

    addpath('Z:\cavity-qed-lvis\Code_calculations\matlab_func\');

    %Remote
    directory = 'Z:\cavity-qed-lvis\Measure_analysis\quantum_feedback\2012-04-26\';
    
    %Local
    %directory = 'D:\dpc-230-data\2012-07-09\';

else

    directory = '~/Measure_analysis/quantum_feedback/2014-04-26/';

end

%Sometimes the data directory and the remote directory are the same
dirremote = directory;

%Larmor precession frequency (in MHz)
%frequency data location - concatenation of file address and name
freq_dat_loc = strcat(dirremote, 'center_freq.txt');
freq_start_row = 1;
freq_start_col = 1;
freq_dat = dlmread(freq_dat_loc, '', freq_start_row, freq_start_col);

freq = freq_dat(:, 1);

freq_err = freq_dat(:, 2);

datafiles=[1];

for dataset = datafiles
       
    data = readg2(directory, dataset);
    
    if (findpulse == 1)
        g2max = max(data(:, 3));
        n = 1;
        for w = 1 : length(data)
            if (data(w, 3) < 0.3 * g2max)
            pulse(n, 1) = data(w, 1);
            pulse(n, 2) = data(w, 3);
            n = n + 1;
            end
        end;
        
        clear g2max;
        
        clear n;
        
        tini = pulse(end, 1) + delay;
        
        pulsewidth = pulse(end, 1) - pulse(1, 1);
        
        disp('Pulse Width: ');
        disp(pulsewidth);
    end
    
    %Isolate the left half of the data
    datasize = size(data, 1);    

    t = data(datasize, 1) - data(1, 1);

    data_reversed = flipud(data);
    
    leftg2 = data_reversed(floor(datasize / 2) : end, 3);

    leftg2_err = data_reversed(floor(datasize / 2) : end, 4);
    
    Llg2 = size(leftg2, 1);

    %Isolate the right half of the data
    rightg2 = data(floor(datasize / 2) : end, 3);

    rightg2_err = data(floor(datasize / 2) : end, 4);
    
    Lrg2 = size(rightg2, 1);
        
    %Select normal and revived curves
    normal = leftg2(floor(Llg2 * (2 * tini / t)) : Llg2);
    
    normal_err = leftg2_err(floor(Llg2 * (2 * tini / t)) : Llg2);

    revived = rightg2(floor(Lrg2 * (2 * tini /t)) : Lrg2);
    
    revived_err = rightg2_err(floor(Lrg2 * (2 * tini /t)) : Lrg2);
    
    size_rightg2 = size(revived, 1);

    start_frac = floor(Lrg2 * (2 * tini / t));
    
    % Approximate "distance" between peaks in number of bins
    min_peak_dist = ceil(min_peak_frac * datasize * (1 / freq(dataset)) / t);
    
    %Center the data about zero
    if (option == 0)
        
        yoff_ini = 0;
        
        yoff_ini_err = 0;
        
        rel_curve = 0;
        
        rel_curve_err = 0;
        
    else
        
        yoff_ini = 1.0;
        
        yoff_ini_err = 0;
        
        if (option == 1)
            
            rel_curve = 0;
            
            rel_curve_err = 0;
            
        else if (option == 2)
                
                [revived_off, revived_off_err] = centerdata(revived, revived_err, min_peak_dist);
                
                [normal_off, normal_off_err] = centerdata(normal, normal_err, min_peak_dist);
                
                rel_curve = revived_off - normal_off;
                
                rel_curve_err = sqrt(revived_off_err.^2 + normal_off_err.^2);
            end
            
        end
        
    end
       
    normal = leftg2(floor(Llg2 * (2 * tini / t)) : floor(Llg2 * ( 2 * tfin / t)));
    
    normal_err = leftg2_err(floor(Llg2 * (2 * tini / t)) : floor(Llg2 * ( 2 * tfin / t)));

    Lnormal = size(normal, 1);
    
    rel_curve = rel_curve(1 : Lnormal);
        
    rel_curve_err = rel_curve_err(1 : Lnormal);
    
    %Compute the differences between leftg2 and all possible translations and
    %scales of rightg2
    unscaleddiff = sum((revived(1 : Lnormal) - normal).^2 ./ (normal_err.^2 + revived_err(1 : Lnormal).^2));

    %Error estimation for phase shift
    
    %find peaks and valleys
    [norm_top, norm_pos_top] = findpeaks(normal, 'MINPEAKDISTANCE', min_peak_dist);

    [norm_bot, norm_pos_bot] = findpeaks(-normal, 'MINPEAKDISTANCE', min_peak_dist);
    
    %concatenate all positions and sort them in ascending order
    normal_peak_pos = sort(cat(1, norm_pos_bot, norm_pos_top));
    
    %Last peak seemed problematic sometimes, so I removed it. This
    %shouldn't affect the outcome negatively.
    normal_peak_pos = normal_peak_pos(1 : end - 1);
    
    %Number of peaks
    normal_num_peaks = size(normal_peak_pos, 1);
    
    %normal inverse time error
    norm_inv_terr = 0;
    
    for inn = 1 : normal_num_peaks - 1
          
        %calculates integer range in between each peak/valley position
        norm_pos_range = normal_peak_pos(inn) : normal_peak_pos(inn + 1);
        
        if (normal_peak_pos(inn + 1) - normal_peak_pos(inn) <= 1)
            continue;
        end
        
        %Fits a line to each peak-valley range
        [norm_par, norm_par_err] = polyfitweighted(norm_pos_range', normal(norm_pos_range), normal_err(norm_pos_range), 1);
        [norm_lin_fit, norm_lin_fit_err] = polyvalweighted(norm_par, norm_par_err, norm_pos_range);

        %Error at the center of the range
        norm_center_pos = floor((normal_peak_pos(inn + 1) - normal_peak_pos(inn)) / 2);
        norm_ctr_terr = abs(min_time_resolve * normal_err(norm_center_pos) / norm_par(1)) / 2.5;
        
        norm_inv_terr = norm_inv_terr + 1 / norm_ctr_terr.^2;
        
        %plot straight lines for error estimation (just for debugging)
        plot(norm_pos_range, norm_lin_fit, 'Color', 'green');
        
    end
    
    final_norm_terr = sqrt(1 / norm_inv_terr);
    
    %find peaks and valleys
    [rev_top, rev_pos_top] = findpeaks(revived(1 : Lnormal), 'MINPEAKDISTANCE', min_peak_dist);

    [rev_bot, rev_pos_bot] = findpeaks(-revived(1 : Lnormal), 'MINPEAKDISTANCE', min_peak_dist);
    
    %concatenate all positions and sort them in ascending order
    revived_peak_pos = sort(cat(1, rev_pos_bot, rev_pos_top));
    
    %Number of peaks
    revived_num_peaks = size(revived_peak_pos, 1);
    
    %revived inverse time error
    rev_inv_terr = 0;
    
    for inr = 1 : revived_num_peaks - 1
        
        %calculates integer range in between each peak/valley position
        rev_pos_range = revived_peak_pos(inr) : revived_peak_pos(inr + 1);
        
        if (revived_peak_pos(inr + 1) - revived_peak_pos(inr) <= 1)
            continue;
        end
        
        %Fits a line to each peak-valley range
        [rev_par, rev_par_err] = polyfitweighted(rev_pos_range', revived(rev_pos_range), revived_err(rev_pos_range), 1);
        [rev_lin_fit, rev_lin_fit_err] = polyvalweighted(rev_par, rev_par_err, rev_pos_range);  
        
        %Error at the center of the range
        rev_center_pos = floor((revived_peak_pos(inr + 1) - revived_peak_pos(inr)) / 2);
        %This error is propagated from the amplitude error of the g2
        %function by a geometrical argument that takes the error
        %bar as a valid distance. There's a rule of thumb that states that
        %to get a better estimate for the error, based on geometrical
        %arguments one has to divide by 5. Here we divide by 2.5, because
        %we only took half the error bar.
        rev_ctr_terr = abs(min_time_resolve * revived_err(rev_center_pos) / rev_par(1)) / 2.5;
        
        rev_inv_terr = rev_inv_terr + 1 / rev_ctr_terr.^2;
        
        %plot straight lines for error estimation (just for debugging)
        plot(rev_pos_range, rev_lin_fit, 'Color', 'red');
        
    end
    
    final_rev_terr = sqrt(1 / rev_inv_terr);

    final_terr = sqrt(final_rev_terr^2 + final_norm_terr^2);
    
    shift_index = 1;
    
    trans_range = floor(mean(diff(rev_pos_top)));
    
    for tshift = 0 : trans_range
        
        rightg2trans = rightg2(start_frac + tshift : size_rightg2 + tshift);
        
        rightg2trans_err = rightg2_err(start_frac + tshift : size_rightg2 + tshift);

        new_size_rightg2 = size(rightg2trans, 1);
        
        trimmed = rightg2trans(1 : Lnormal);
        
        trimmed_err = rightg2trans_err(1 : Lnormal);
     
        for j = 1 : ny;
                
            yscaled = (trimmed - yoff_ini - rel_curve) * scy(j) + yoff_ini;
            
            yscaled_err = sqrt((trimmed_err * scy(j)).^2 + (rel_curve_err * scy(j)).^2 + (yoff_ini_err * (1 - scy(j)))^2);
                
            ydiff(j, [1 2]) = [sum((yscaled - normal).^2 ./ (yscaled_err.^2 + normal_err.^2)) scy(j)];

        end

        %search diff for the optimal value of scy
        [ymin, ymin_index] = min(ydiff(:, 1));

           
        %I don't think this index on scy_opt is necessary...revisit later
        %Edit: More pressing tasks have prevented me from revisiting this
        %question, but it is just a question of efficiency, not accuracy.
        scy_opt(shift_index) = ydiff(ymin_index, 2);
        
        yscaled_min = scy_opt(shift_index) * (trimmed - yoff_ini - rel_curve) + yoff_ini;
        
        yscaled_min_err = sqrt((trimmed_err * scy_opt(shift_index)).^2 + (rel_curve_err * scy_opt(shift_index)).^2 + (yoff_ini_err * (1 - scy_opt(shift_index)))^2);
            
        xydiff(shift_index, [1 2 3]) = [sum((yscaled_min - normal).^2 ./ (yscaled_min_err.^2 + normal_err.^2)) tshift scy_opt(shift_index)];

        shift_index = shift_index + 1;
        
    end

    %Search for the optimal translation and scy
    [xymin, xymin_index] = min(xydiff(:, 1));

    shift_best = xydiff(xymin_index, 2);
    
    phase_best = shift_best * min_time_resolve * freq(dataset);
    
    %Propagate the errors in time difference to phase difference
    phase_err = sqrt((freq(dataset) * final_terr).^2 + (freq_err(dataset) * shift_best * min_time_resolve).^2);
    
    scy_best = xydiff(xymin_index, 3);
    
    [rightg2_opt, rightg2_opt_err] = g2_transform(rightg2, Lnormal, size_rightg2, rightg2_err, start_frac, [shift_best scy_best], yoff_ini, yoff_ini_err, rel_curve, rel_curve_err);
    
    disp(strcat('data # ', int2str(dataset), ' fractional improvement:'));
    
    disp(xymin / unscaleddiff);
    
    disp('Reduced Chi2');
    
    disp(xymin / (size(rightg2_opt, 1) - fit_param));
    
    disp('Best parameters:');
    
    disp(shift_best);
    
    disp(scy_best);
    
    disp('Best phase difference:');
    
    disp(phase_best);
    
    %Checking error in y scaling factor
    frac_var = 0.11;
    
    num_elements = 5000;
    
    tolerance = 1e-5;
    
    %Define variation range around best value
    scy_var = linspace(scy_best * (1 - frac_var), scy_best * (1 + frac_var), num_elements);
    
    %Preallocation for speed
    chi2_var_scy = zeros(1, num_elements);
    scy_possible = zeros(1, 1);
    
    for kk = 1 : num_elements
    
        chi2_var_scy(kk) = chi2diff(normal, rightg2, Lnormal, size_rightg2, normal_err, rightg2_err, start_frac, [shift_best scy_var(kk)], yoff_ini, yoff_ini_err, rel_curve, rel_curve_err);
        
    end
    
    polyelements = 1000000;
    
    poly_scy = polyfit(scy_var, chi2_var_scy, 2);
    
    scy_range = linspace(scy_best * (1 - frac_var), scy_best * (1 + frac_var), polyelements);
    
    chi2_poly_scy = polyval(poly_scy, scy_range);
    
    for kk = 1 : polyelements
     
         if (chi2_poly_scy(kk) > xymin + 1 - tolerance) && (chi2_poly_scy(kk) < xymin + 1 + tolerance)
            
             scy_possible(end + 1) =  scy_range(kk) - scy_best;
            
         end
         
     end
    
    %Vectorized alternative to the loop commented out above
        %mask = chi2_poly_scy > xymin + 1;
        %[chi, pos] = min(chi2_poly_scy(mask));
        %scy_possible = scy_range(pos) - scy_best;
    
    %Test of best scy
    [min_fit, scy_min_fit_index] = min(chi2_poly_scy);
    
    disp('discrepancy between polyfit and previous minimization algorithm: ');
    
    disp((scy_range(scy_min_fit_index) - scy_best) / scy_best);
    
    disp('Errors yscaling: ');
    
    disp(scy_possible);
    %disp(max(scy_possible));
    %disp(min(scy_possible));
    
    disp('phase error estimate:');
    
    disp(phase_err);
    
    
    
    scaled_plots = plot(rightg2_opt);
    
    %random rgb code for setting color plots 
    rgb_codes = rand(1, 3);
    
    set(scaled_plots, 'Color', rgb_codes);

end

normal_plot = plot(normal);

set(normal_plot, 'Color', 'black');

figure(gcf); %bring plot to front

plot(revived(1 : Lnormal));

[tops, ptops] = findpeaks(revived(1 : Lnormal), 'MINPEAKDISTANCE', ceil(datasize / 100));

%If the first entry is a peak, then include it
if (max(revived(1 : Lnormal)) == revived(1))
    for w = length(tops) : 1
        tops(w + 1) = tops(w);
        ptops(w + 1) = ptops(w);
    end
    tops(1) = revived(1);
    ptops(1) = 1; 
end

[bots, pbots] = findpeaks(-revived(1 : Lnormal), 'MINPEAKDISTANCE', ceil(datasize / 100));

bots = -bots;

%If the first entry is a peak, then include it
if (min(revived(1 : Lnormal)) == revived(1))
    for w = length(tops) : 1
        bots(w + 1) = bots(w);
        pbots(w + 1) = pbots(w);
    end
    bots(1) = revived(1);
    pbots(1) = 1; 
end

top_coeff = polyfit(ptops, tops, min(2, length(ptops) - 1));
    
topfit = polyval(top_coeff, 1 : Lnormal);
    
bot_coeff = polyfit(pbots, bots, min(2, length(pbots) - 1));
    
botfit = polyval(bot_coeff, 1 : Lnormal);

plot(topfit);

plot(botfit);

hold off

figure();

hold on

plot(scy_var, chi2_var_scy); 

plot(scy_range, chi2_poly_scy, 'Color', 'red');

hold off