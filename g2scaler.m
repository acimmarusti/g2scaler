clear all
hold off

%range of the t<0 data to be analyzed (in microseconds)
tini = 2.7;
tfin = 4;

%how many different scales to check
n = 100;

%translation range (CAREFUL!)
trans_range = 100;
half_range = floor(trans_range / 2);

%range of scales to check
scx = linspace(0.5, 2, n);
scy = linspace(1, 0, n);

%real data or artificial data (for testing)?
realdata = true;

datafiles=[12 12];

hold on

for dataset = datafiles(1) : datafiles(2)
    
    if realdata
        
    %Windows path
    directory = 'Z:\cavity-qed-lvis\Measure_analysis\quantum_feedback\2012-03-23\';

    %GNU/Linux path
    %directory = '/home/candres/Measure_analysis/quantum_feedback/2012-03-23/';

    %data location - concatenation of file address and name
    datlocation = strcat(directory, 'data', int2str(dataset),'_g2.asc');
    
    %loads data from start_row/start_col and beyond (whitespace delimeter)
    start_row = 9;
    start_col = 0;
    data = dlmread(datlocation, '', start_row, start_col);
    
    else
        x = linspace(-16, 16, 2000);
        Lx = size(x, 2);
        data = [x; diag(zeros(2000))'; - sin(3 * x(1 : ceil(Lx / 2))), 0.8 * sin((2.5 * x(ceil(Lx / 2) + 1 : end) + 0.5))]';
    end
    
    %Average error per dataset
    error = mean(data(:, 4));
    
    %Isolate the left half of the data
    datasize = size(data, 1);
    
    t = data(datasize, 1) - data(1, 1);

    leftg2 = flipud(data);

    leftg2 = leftg2(floor(datasize / 2) : end, 3);

    Llg2 = size(leftg2, 1);

    %Isolate the right half of the data
    rightg2 = data(floor(datasize / 2) : end, 3);

    Lrg2 = size(rightg2, 1);

    normal = leftg2(floor(Llg2 * (2 * tini / t)) : floor(Llg2 * ( 2 * tfin / t)));

    Lnormal = size(normal, 1);
    
    revived = rightg2(floor(Lrg2 * (2 * tini / t)) : Lrg2);
    
    size_rightg2 = size(revived, 1);

    start_frac = floor(Lrg2 * (2 * tini / t));
    
    %Compute the differences between leftg2 and all possible translations and
    %scales of rightg2
    initial_param = [0 1 1];
    
    unscaleddiff = sum(dataset_diff(normal, rightg2, Lnormal, size_rightg2, error, error, start_frac, initial_param).^2);
    
    %shift xscale yscale. The function must accept this vector as input
    %The function must only output a vector of the differences. The
    %minimization Routine will square and sum them.
    guess_param = [3 1 0.95];
    
    upper_bounds = [5 1.5 1.5];
    lower_bounds = [-5 0.25 0.25];
    
    %The optimization routine will return a vector with the optimal
    %parameters and optionally the minimum chi^2
    [best_param, min_chi2] = lsqnonlin(@(x)dataset_diff(normal, rightg2, Lnormal, size_rightg2, error, error, start_frac, x), guess_param, lower_bounds, upper_bounds);
    
    disp(strcat('data # ', int2str(dataset), ' fractional improvement:'));
    
    disp(min_chi2 / unscaleddiff);
    
    disp(best_param);
    
    rightg2_opt = rightg2(start_frac + best_param(1) : size_rightg2 + best_param(1));       %shift
    
    new_size_rightg2 = size(rightg2_opt, 1);
    
    rightg2_opt = best_param(3) * interpft(rightg2_opt, ceil(new_size_rightg2 * best_param(2)));% + yoffset_best;  	 %scale
    
    rightg2_opt = rightg2_opt(1 : Lnormal);                                             		 %trim
    
    scaled_plots = plot(rightg2_opt);
    
    %random rgb code for setting color plots 
    rgb_codes = rand(1, 3);
    
    set(scaled_plots, 'Color', rgb_codes);

end

normal_plot = plot(normal);

set(normal_plot, 'Color', 'black');

figure(gcf); %bring plot to front

%axis([0 Lnormal 0.8 1.4]);

plot(revived(1 : Lnormal));

hold off
