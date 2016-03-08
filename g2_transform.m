% Computes transformations on data2 for comparison with data1
% Data2 can undergo transformations such as:
% - Starting point (start)
% - Horizontal index shift (par1)
% - Vertical scaling (par2)
% Syntax is:
% g2_transform(data1, data2, size_dat1, size_dat2, data1_err, data2_err, start, par)

function [transformed, transformed_err] = g2_transform(data2, size_dat1, size_dat2, data2_err, start, par, offset, offset_err, off_curve, off_curve_err)

translated = data2(floor(start + par(1)) : floor(size_dat2 + par(1)));

translated_err = data2_err(floor(start + par(1)) : floor(size_dat2 + par(1)));

cut = translated(1 : size_dat1);

cut_err = translated_err(1 : size_dat1);

transformed = par(2) * (cut - offset - off_curve) + offset;

transformed_err = sqrt((cut_err * par(2)).^2 + (off_curve_err * par(2)).^2 + (offset_err * (1 - par(2)))^2);

end

