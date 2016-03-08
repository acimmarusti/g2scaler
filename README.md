### Brute-force data-to-data fitting.

g2scaler_brute is a fitting routine designed to extract the phase, frequency, and amplitude difference from two data sets. It is specifically designed to operate on second order correlation data from our Quantum Feedback in Cavity QED experiment.

It achieves this by fixing the scale of one data vector (in this case, the left half of the g2), and transforming the vector until they are optimally matched in a least-squares sense.

The program performs optimization according to the following procedure. It originally also performed horizontal scaling, but early results indicated that a frequency shift was not present.

The program imports the specified data and performs pre-processing. It separates the left and right halves of the data, the left being the "unpulsed" target and the right being the data to transform to fit. It identifies the region of iterest after a pulse and isolates the range in both sets of data. The corresponding errors in the data are stored alongside the target and fit data.

Using a custom function called centerdata, the program removes the best quadratic trend from the transformed data and centers it vertically about zero. This isolates the fitting from stray experimental effects independent of amplitude and phase. The errors are propagated to reflect the operation. This routine uses custom functions called polyfitweighted and polyvalweighted. Polyfitweighted performs error-weighted polynomial fits and returns fitting parameters with propagated errors. Polyvalweighted generates a vector of y values and y errors from the fit object returned by Polyfitweighted.

The syntax is documented in the code. An external loop iterates over all possible horizontal translations by shifting the elements of the transformed vector. Every iteration the sum of squared differences is computed and saved. A nested loop iterates over all possible vertical scales. Every iteration the sum of squared differences is computed and saved.

At each iteration of this brute-force search, errors are propagated and stored to reflect each transformation.

When the looping is complete the program searches the results of every transformation for the smallest difference and stores the corresponding xshift-yscale parameters.

To find the error in the parameters we extracted via brute force, we perform a chi squared analysis individually on each. Using a custom function called chi2diff, we evaluate the chi^2 goodness-of-fit in the neighborhood of the optimal parameter. This yields a parabola that can be easily fit.The sharpness of the parabola corresponds to the parameter error bounds, and the location of the minima corresponds to the optimal value of the parameter. This optimal value is very near what is found via brute force, and serves as both a refinement and a good consistency check.