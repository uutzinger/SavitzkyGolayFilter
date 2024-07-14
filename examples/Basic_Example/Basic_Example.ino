/**
   SavLayFilter Basic Example
   Creates a generic sin wave and adds noise from random numbers
   Two filters are then run on it with a large windowsize and a small windowsize
   It is best to open Serial plotter for best results

   For non-Teensy users, elapsedMillis is required
   It can be downloaded from: https://github.com/pfeerick/elapsedMillis/archive/master.zip
   
   Author: James Deromedi
   Updates: Urs Utzinger
   License: MIT License
**/


#include <Arduino.h>
#include "SavLayFilter.h"

#undef SMALL               // max window size is 25 if SMALL is defined, otherwise 63

// Window size, order, and derivative
#define WINDOW_SIZE_1 5    // must be odd 3...25 or 63 
#define WINDOW_SIZE_2 25   // must be odd 3...25 or 63
#define ORDER 2            // 1 = Linear, 2 = Quadratic, 3 = Cubic, 4 = Quartic, 5 = Quintic
#define DERIVATIVE 0       // 0 = Smooth, 1 = First Derivative, 2 = Second Derivative

float phase = 0.0;
float twopi = 3.14159 * 2;

int32_t signalValue;
unsigned long lastTime;
unsigned long currentTime;

// Initialize the SavLayFilter
SavLayFilter smallFilter(WINDOW_SIZE_1, ORDER, DERIVATIVE);
SavLayFilter largeFilter(WINDOW_SIZE_2, ORDER, DERIVATIVE);

//==================================================================================================

void setup() {
  // Start the serial communication
  Serial.begin(192000);
  
  // Wait for serial port to connect
  while (!Serial) {
    ; // wait for serial port to connect. Needed for native USB
  }

  lastTime = micros();
}

//==================================================================================================

void loop() {

  currentTime = micros();

  if ((currentTime - lastTime) > 3000){
    lastTime = currentTime;
    // Generate a random value
    int32_t randomValue = random(0, 500);

    signalValue = int32_t(sin(phase) * 1000.0 + 2000.0) + randomValue;   // creates sin wave pattern with A = 1000 and shifted up by 2000
    phase = phase + 0.02;                                     // propaget the sine wave
    if (phase >= twopi) phase = 0;                            // resets the phase

    // Update the filter with the random value
    int32_t filteredValue_small = smallFilter.update(signalValue);
    int32_t filteredValue_large = largeFilter.update(signalValue);

    Serial.print(signalValue);                                //Raw Value [Blue line]
    Serial.print(",");
    Serial.print(filteredValue_small);                      //Smoothed value of smaller window [Orange line]
    Serial.print(",");
    Serial.println(filteredValue_large);                    //Smoothed value of smaller window [Red line]

  } else {
    delay(0);
  }

}
