/**
   SavitzkyGolayFilter Derivative Example
   Creates a generic sin wave
   A filter is then run on it to compute the first derivative or slope of the sin wave
   The derivative is normalized by multipling by the amplitude of the sin wave
   There is a known issue that when the sin wave completes a period there is a dip in the derivative
   It is best to open Seral Plotter for best results

   Author: James Deromedi
   Updates: Urs Utzinger
   License: MIT License
**/

#include <Arduino.h>
#include "SavitzkyGolayFilter.h"

#define SMALL             // max window size is 25, otherwise 63

// Window size, order, and derivative
#define WINDOW_SIZE 5    // must be odd 3...25 or 63 
#define ORDER 2          // 1 = Linear, 2 = Quadratic, 3 = Cubic, 4 = Quartic, 5 = Quintic
#define DERIVATIVE 1     // 0 = Smooth, 1 = First Derivative, 2 = Second Derivative

float phase = 0.0;
float phase_inc = 0.02;
float phase_inc_factor = 1.5;
float twopi = 3.14159 * 2;

int32_t signalValue;
unsigned long lastTime;
unsigned long currentTime;

// Initialize the SavLayFilter
SavLayFilter firstDerivative(WINDOW_SIZE, ORDER, DERIVATIVE);

void setup() {
  // Start the serial communication
  Serial.begin(500000);
  delay(2000);
  Serial.println("Starting SavGol");
  lastTime = micros();
}

void loop() {
  currentTime = micros();

  if (currentTime-lastTime > 3000){
    lastTime = currentTime;
    // Generate a random value
    int32_t randomValue = random(0, 500);

    signalValue = int32_t(sin(phase) * 1000.0 + 2000.0) + randomValue;   // creates sin wave pattern with A = 1000 and shifted up by 2000
    phase = phase + phase_inc;                                     // propaget the sine wave
    if (phase >= twopi) {
      phase = 0;                            // resets the phase
      phase_inc = phase_inc * phase_inc_factor;
      if (phase_inc > twopi/5) {
        phase_inc_factor = 0.66;
        phase_inc = twopi/5.;
      } else if (phase_inc < 0.02) {
        phase_inc_factor = 1.5;
        phase_inc = 0.02;
      }
    }

    // Update the filter with the random value
    int32_t filteredValue = firstDerivative.update(signalValue);

    Serial.print(signalValue);                                //Raw Value [Blue line]
    Serial.print(",");
    Serial.println(filteredValue*10);                         //Smoothed value of smaller window [Orange line]

  } else {
    delay(0);
  }

}
