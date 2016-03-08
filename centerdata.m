% Computes the 'average' (center) curve from the envelope
% Data must be oscillatory
% Syntax is:
% centerdata(data, data_err, min_dist_peak)
function [center, center_err] = centerdata(data, data_err, min_dist_peak)

    ldata = length(data);
        
    [tops, ptops] = findpeaks(data(1 : ldata), 'MINPEAKDISTANCE', min_dist_peak);
    
    %If the first entry is a peak, then include it
    if (max(data) == data(1))
        for w = length(tops) : 1
            tops(w + 1) = tops(w);
            ptops(w + 1) = ptops(w);
        end
        tops(1) = data(1);
        ptops(1) = 1; 
    end
        
    [bots, pbots] = findpeaks(-data(1 : ldata), 'MINPEAKDISTANCE', min_dist_peak);
    bots = - bots;

    %If the first entry is a peak, then include it
    if (min(data) == data(1))
        for w = length(tops) : 1
            bots(w + 1) = bots(w);
            pbots(w + 1) = pbots(w);
        end
        bots(1) = data(1);
        pbots(1) = 1; 
    end
    
    tops_err = data_err(ptops);
    
    bots_err = data_err(pbots);
    
    top_poly_order = min(2, length(ptops) - 1);

    bot_poly_order = min(2, length(pbots) - 1);
    
    %Weighted polynomial fits and evaluation
    [top_fit_par, top_fit_par_err] = polyfitweighted(ptops, tops, tops_err, top_poly_order);
    
    [bot_fit_par, bot_fit_par_err] = polyfitweighted(pbots, bots, bots_err, bot_poly_order);
    
    [topfit, topfit_err] = polyvalweighted(top_fit_par, top_fit_par_err, 1, ldata);
    
    [botfit, botfit_err] = polyvalweighted(bot_fit_par, bot_fit_par_err, 1, ldata);

    %Average to find center curve
    center = (topfit + botfit)' / 2;
    
    center_err = (sqrt(topfit_err.^2 + botfit_err.^2))' / 2;
end