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
binsize = 10;

%Minumum time resolution (in microseconds)
min_time_res = 164.15e-6;

%Larmor precession frequency (in MHz)
larmor_freq = 2.5;

%total fitting parameters:
%yscaling, xshift, upper envelope(3) and lower envelope(3)
fit_param = 8;

datafiles=[11 11];

min_time_resolve = binsize * min_time_res;

hold on

for dataset = datafiles(1) : datafiles(2)
        
    %Windows path remote folder
    %directory = 'Z:\cavity-qed-lvis\Measure_analysis\quantum_feedback\2012-04-26\';
    
    %Windows path local folder
    directory = 'D:\dpc-230-data\2012-05-16\';

    %GNU/Linux path
    %directory = '/home/candres/Measure_analysis/quantum_feedback/2012-04-26/';

    %data location - concatenation of file address and name
    datlocation = strcat(directory, 'data', int2str(dataset),'_g2.asc');
    
    %loads data from start_row/start_col and beyond (whitespace delimeter)
    start_row = 9;
    start_col = 0;
    data = dlmread(datlocation, '', start_row, start_col);
    
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
                
                [revived_off, revived_off_err] = centerdata(revived, revived_err, datasize);
                
                [normal_off, normal_off_err] = centerdata(normal, normal_err, datasize);
                
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
    
    [norm_top, norm_pos_top] = findpeaks(normal, 'MINPEAKDISTANCE', ceil(datasize/100));

    [norm_bot, norm_pos_bot] = findpeaks(-normal, 'MINPEAKDISTANCE', ceil(datasize/100));
    
    normal_peak_pos = sort(cat(1, norm_pos_bot, norm_pos_top));
    
    %Last peak seemed problematic so I removed it
    normal_peak_pos = normal_peak_pos(1 : end - 1);
    
    normal_num_peaks = size(normal_peak_pos, 1);
    
    norm_inv_terr = 0;
    
    for inn = 1 : normal_num_peaks - 1
        
        norm_pos_range = normal_peak_pos(inn) : normal_peak_pos(inn + 1);
        
        norm_par = polyfit(norm_pos_range', normal(norm_pos_range), 1);

        %Error at the center of the range
        norm_center_pos = floor((normal_peak_pos(inn + 1) - normal_peak_pos(inn)) / 2);
        norm_ctr_terr = abs(min_time_resolve * normal_err(norm_center_pos) / norm_par(1));
        
        norm_inv_terr = norm_inv_terr + 1 / norm_ctr_terr.^2;
        
        %plot straight lines for error estimation (just for debugging)
        norm_lin_fit = polyval(norm_par, norm_pos_range);
        plot(norm_pos_range, norm_lin_fit, 'Color', 'green');
        
    end
    
    final_norm_terr = sqrt(1 / norm_inv_terr);
    
    [rev_top, rev_pos_top] = findpeaks(revived(1 : Lnormal), 'MINPEAKDISTANCE', ceil(datasize/100));

    [rev_bot, rev_pos_bot] = findpeaks(-revived(1 : Lnormal), 'MINPEAKDISTANCE', ceil(datasize/100));
    
    revived_peak_pos = sort(cat(1, rev_pos_bot, rev_pos_top));
    
    revived_num_peaks = size(revived_peak_pos, 1);
    
    rev_inv_terr = 0;
    
    for inr = 1 : revived_num_peaks - 1
        
        rev_pos_range = revived_peak_pos(inr) : revived_peak_pos(inr + 1);
        
        rev_par = polyfit(rev_pos_range', revived(rev_pos_range), 1);
        
        %Error at the center of the range
        rev_center_pos = floor((revived_peak_pos(inr + 1) - revived_peak_pos(inr)) / 2);
        rev_ctr_terr = abs(min_time_resolve * revived_err(rev_center_pos) / rev_par(1));
        
        rev_inv_terr = rev_inv_terr + 1 / rev_ctr_terr.^2;
        
        %plot straight lines for error estimation (just for debugging)
        rev_lin_fit = polyval(rev_par, rev_pos_range);  
        plot(rev_pos_range, rev_lin_fit, 'Color', 'red');
        
    end
    
    final_rev_terr = sqrt(1 / rev_inv_terr);

    %Propagate the errors in time difference to phase difference
    phase_err = larmor_freq * sqrt(final_rev_terr^2 + final_norm_terr^2);
    
    
    
    
        trimmed = rightg2(start_frac : Lnormal);
        
        trimmed_err = rightg2_err(start_frac : Lnormal);
     
        for j = 1 : ny;
                
            yscaled = (trimmed - yoff_ini - rel_curve) * scy(j) + yoff_ini;
            
            yscaled_err = sqrt((trimmed_err * scy(j)).^2 + (rel_curve_err * scy(j)).^2 + (yoff_ini_err * (1 - scy(j)))^2);
                
            ydiff(j, [1 2]) = [sum((yscaled - normal).^2 ./ (yscaled_err.^2 + normal_err.^2)) scy(j)];

        end

        %search diff for the optimal value of scy
        [ymin, ymin_index] = min(ydiff(:, 1));


    %Search for the optimal translation and scy
    [xymin, xymin_index] = min(xydiff(:, 1));

    %phase_best = shift_best * min_time_resolve * larmor_freq;
    
    scy_best = ymin;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I'm leaving so I do not yet know if the below code works with the changes I have made%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
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
    frac_var = 0.1;
    
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

% [tops, ptops] = findpeaks( revived(1 : end), 'MINPEAKDISTANCE', ceil(datasize/100));
%     %If the first entry is a peak, then include it
%     if(max(revived(1:Lnormal))==revived(1))
%         for w=length(tops):1
%             tops(w+1)  =  tops(w);
%             ptops(w+1) = ptops(w);
%         end
%         tops(1) = revived(1);
%         ptops(1)= 1; 
%     end
% [bots, pbots] = findpeaks(-revived(1 : end), 'MINPEAKDISTANCE', ceil(datasize/100));
% bots = -bots;
%     %If the first entry is a peak, then include it
%     if(min(revived(1:Lnormal))==revived(1))
%         for w=length(tops):1
%             bots(w+1)  =  bots(w);
%             pbots(w+1) = pbots(w);
%         end
%         bots(1) = revived(1);
%         pbots(1)= 1; 
%     end
% 
% %Quadratic fit (linear if less than three peaks are detected)
% topfit = polyval(polyfit(ptops, tops, min(2,length(ptops)-1)), 0 : length(revived));
% botfit = polyval(polyfit(pbots, bots, min(2,length(pbots)-1)), 0 : length(revived));
% 
% 
% %Since we subtract the envelope in the beginning, it doesn't
% %make sense to re-calculate and plot it here. We should find a way to plot
% %it the first time it is calculated, however
% %plot(topfit);
% %plot(botfit); 

plot(revived(1 : Lnormal));

hold off

figure();

hold on

plot(scy_var, chi2_var_scy); 

plot(scy_range, chi2_poly_scy, 'Color', 'red');

hold off