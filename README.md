# SavitzkyGolayFilter
This library contains provides one dimensional Savitzky-Golay's filtering algorithm

Documentation can be found within the code extras folder.

The orginal code was developed by James Deromedi.

This fork was created by Urs Utzinger.

It provides python code to compute the filter coefficients, so you can extend the filter kernels.
SciPy signal module was used. Most recreated kernels match the original published ones. 
However some do not seem to be calculated correctly in SciPy, e.g.: quartic first derivative, quintic and sexic first derivative, quintic second derivative.
Some of the scaling factors in the original paper are very large and not useful for integer math e.g. 4E8.

It provide a small and large version of the filter tables for shorter and longer window sizes. This results in smaller or larger memory use.

The convolution was improved for speed and simplicity.
A ring buffer was implemented and data is assumed to be int32 instead of double. The kernel is unfolded so that convolution is simpler and faster.
