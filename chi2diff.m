% Function that computes Chi^2 between two datasets (arrays)
% Data2 can undergo transformations such as:
% - Starting point (start)
% - Horizontal index shift (par1)
% - Vertical scaling (par2)
% Syntax is:
% chi2diff(data1, data2, size_dat1, size_dat2, data1_err, data2_err, start, par, offset, offset_err, off_curve, off_curve_err)

function output = chi2diff(data1, data2, size_dat1, size_dat2, data1_err, data2_err, start, par, offset, offset_err, off_curve, off_curve_err)

differences = dataset_diff(data1, data2, size_dat1, size_dat2, data1_err, data2_err, start, par, offset, offset_err, off_curve, off_curve_err);

output = sum(differences.^2);