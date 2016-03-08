% Function that computes differences between two datasets (arrays)
% Data2 can undergo transformations such as:
% - Starting point (start)
% - Horizontal index shift (par1)
% - Vertical scaling (par2)
% Syntax is:
% dataset_diff(data1, data2, size_dat1, size_dat2, data1_err, data2_err, start, par, offset, offset_err, off_curve, off_curve_err)

function differences = dataset_diff(data1, data2, size_dat1, size_dat2, data1_err, data2_err, start, par, offset, offset_err, off_curve, off_curve_err)

[shrunk, shrunk_err] = g2_transform(data2, size_dat1, size_dat2, data2_err, start, par, offset, offset_err, off_curve, off_curve_err);

differences = (shrunk - data1) ./ sqrt(data1_err.^2 + shrunk_err.^2);