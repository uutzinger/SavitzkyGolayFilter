# SavitzkyGolayFilter

This library provides one dimensional Savitzky-Golay's filtering algorithm.
This filter better smooths data with peaks and valleys such as ECG or spectra.

Example smoothing and differentiation is shown below.

<img src="assets/SavGolSmooth.png" alt="Smooth" width="800"/>

Diagram showing short and long windowlength using 2nd order polynomia smoothing.
There is signal delay matching window length and decrease of fidelity at high frequencies.

<img src="assets/SavGolDerivative.png" alt="Derivative" width="800"/>

Diagram showing differentiation using short window length. Derivative increases when osciallation increases. Short window length picks up noise.

## Introduction

SavGol filter approximates a polynomia to the data. A window is shifted over the data and the center point of the window is replaced with the point of the fit.
One can calculate the derivative of the polynomium and compute derivatives instead of smoothed data.

This can be implemented as a convolution, where the kernel is computed in advance.

## Documentation

Documentation can be found within the code extras folder.

## Implementation, Changes

The orginal Arduino code was developed by James Deromedi.
This fork was created by Urs Utzinger.

It provides python code to compute the filter coefficients, so you can extend the library with your own filter kernels.
The SciPy signal module is used. Most recreated kernels match the original published ones. 
The original Savitzky Golay paper had scaling facators up to 4E8. Here they are limitted to a max of 16 bit integer * max window width.

The code provides a small and large version of the filter tables for shorter and longer window sizes. This results in smaller or larger memory use.

The convolution was rewritten for speed and simplicity.
A ring buffer was implemented and data can be int32, float or double. 

One filter update on 25 values takes about 3-4 microseconds on ESP32.

Urs Utzinger,
July 2024
