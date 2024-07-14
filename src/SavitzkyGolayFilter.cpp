/**
 * Arduino Savitzky Golay Library - Version 1.0.1
 * by James Deromedi <jmderomedi@gmail.com>
 *
 * Arduino Savitzky Golay Library - Version 1.3.0
 * Urs Utzinger
 * 
 * This Library is licensed under the MIT License
 */

#if ARDUINO >= 100
  #include "Arduino.h"
#else
  #include "WProgram.h"
#endif

#define SMALL // use the smaller sized tables
#include "SavitzkyGolayFilter.h"
#include <iostream>
#include <vector>

SavLayFilter::SavLayFilter(int windowSize, int order, int derivative):
  _buffer(windowSize), _derivative(derivative), 
  _bufferSize(windowSize), _windowSize(windowSize),
  _order(order),
  _head(0), _isBufferFull(false) 
{
  if (_windowSize > MAX_WINDOW_SIZE) { 
    _windowSize = MAX_WINDOW_SIZE; 
    _buffer.resize(_windowSize); // Resize buffer to MAX_WINDOW_SIZE
  }

  _halfWindowSize = windowSize / 2;
  _kernelPointer = (windowSize - 3) / 2;
  #ifdef SMALL
    if (_kernelPointer > 11) { _kernelPointer = 11; }
  #else
    if (_kernelPointer > 30) { _kernelPointer = 30; }
  #endif
  initializeConvolutionTable(order, derivative);
  _norm = _convolutionTable[_kernelPointer][0];

  // Create the full mirrored kernel
  _mirroredKernel.resize(_windowSize);
  for (int i = 0; i < _halfWindowSize; i++) {
    _mirroredKernel[i] = _convolutionTable[_kernelPointer][i + 1];
    if (_derivative % 2 == 1) {
      _mirroredKernel[_windowSize - 1 - i] = -_convolutionTable[_kernelPointer][i + 1];
    } else {
      _mirroredKernel[_windowSize - 1 - i] = _convolutionTable[_kernelPointer][i + 1];
    }
  }
  _mirroredKernel[_halfWindowSize] = _convolutionTable[_kernelPointer][_halfWindowSize + 1];
}

void SavLayFilter::initializeConvolutionTable(int order, int derivative) {
  if (order == 1) {
    if (derivative == 0) {
      _convolutionTable = _linearSmooth;
    } else if (derivative == 1) {
      _convolutionTable = _linearFirstDerivative;
    }
  } else if (order == 2) {
    if (derivative == 0) {
      _convolutionTable = _quadraticSmooth;
    } else if (derivative == 1) {
      _convolutionTable = _quadraticFirstDerivative;
    } else if (derivative == 2) {
      _convolutionTable = _quadraticSecondDerivative;
    }
  } else if (order == 3) {
    if (derivative == 0) {
      _convolutionTable = _cubicSmooth;
    } else if (derivative == 1) {
      _convolutionTable = _cubicFirstDerivative;
    } else if (derivative == 2) {
      _convolutionTable = _cubicSecondDerivative;
    }
  } else if (order == 4) {
    if (derivative == 0) {
      _convolutionTable = _quarticSmooth;
    } else if (derivative == 1) {
      _convolutionTable = _quarticFirstDerivative;
    } else if (derivative == 2) {
      _convolutionTable = _quarticSecondDerivative;
    }
  } else if (order == 5) {
    if (derivative == 0) {
      _convolutionTable = _quinticSmooth;
    } else if (derivative == 1) {
      _convolutionTable = _quinticFirstDerivative;
    } else if (derivative == 2) {
      _convolutionTable = _quinticSecondDerivative;
    }
  }
}

int32_t SavLayFilter::update(int32_t newValue) {
  // Circular buffer logic
  _buffer[_head] = newValue;
  _head = (_head + 1) % _bufferSize;

  if (!_isBufferFull) {
    // Return new value until buffer is full
    if (_head == 0) {
      _isBufferFull = true;
    }
    return newValue;
  } else {
    
    // Convolve
    _sum = 0;
    for (int i = 0; i < _windowSize; i++) {
      int _index = (_head + i) % _bufferSize;
      _sum += _buffer[_index] * _mirroredKernel[i];
    }

    // Normalize
    _sum /= _norm;

    return _sum;
  }
}

float SavLayFilter::update(float newValue) {
  // Circular buffer logic
  _buffer_float[_head] = newValue;
  _head = (_head + 1) % _bufferSize;

  if (!_isBufferFull) {
    // Return new value until buffer is full
    if (_head == 0) {
      _isBufferFull = true;
    }
    return newValue;
  } else {
    
    // Convolve
    _sum_float = 0;
    for (int i = 0; i < _windowSize; i++) {
      int _index = (_head + i) % _bufferSize;
      _sum_float += _buffer_float[_index] * (float)_mirroredKernel[i];
    }

    // Normalize
    _sum_float /= (float)_norm;

    return _sum_float;
  }
}

double SavLayFilter::update(double newValue) {
  // Circular buffer logic
  _buffer_double[_head] = newValue;
  _head = (_head + 1) % _bufferSize;

  if (!_isBufferFull) {
    // Return new value until buffer is full
    if (_head == 0) {
      _isBufferFull = true;
    }
    return newValue;
  } else {

    // Convolve
    _sum_double = 0;
    for (int i = 0; i < _windowSize; i++) {
      int _index = (_head + i) % _bufferSize;
      _sum_double += _buffer_double[_index] * (double)_mirroredKernel[i];
    }

    // Normalize
    _sum_double /= (double)_norm;

    return _sum_double;
  }
}
