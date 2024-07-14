################################################################################################3
# Compute Savitzky Golay coefficients for various window sizes and polynomial orders
# Uses scipy.signal.savgol_coeffs to generate the coefficients
# Brute force scaling so that (coefficients/scale) are almost same as original floating point coefficients 
# Urs Utzinger
# 2024
################################################################################################3

import numpy as np
from scipy.signal import savgol_coeffs

# quadCubic:     y=ax3+bx2+cx+d
# quarticQuintuc y=ax5+bx4+cx3+dx2+ex+f

# We have 32bit intger math on ESP32 and ESP8266
# We have 32bit floating point math on ESP32 and ESP8266
# We have 64bit intger math on teensy
# We have 64bit floating point math on teensy
# 
# If we use 32bit integer math our scaling should not exceed half the bit depth assuming the measured values do not exceed half the bit depth also

# max_s = np.iinfo(np.int16).max
max_s = 1024000 # Sacitzky-Golay have in their paper scaling coeefccients up to 430 Million, 
                # should not go above 1 Million because of memory restrictions on 16Gb machine

def findScale(coeffs, max_scale=256000):
    # vectorized, will calcualte all residuals for all scaling factors at once
    # limit max scale to about 1Million with 16Gb of memory
    best_scaling = None
    eps = np.finfo(np.float64).eps
    eps10 = 10*eps
    under_thresh_limit = 1.0-eps10

    # vectorization
    scaling_factors = np.arange(1, max_scale)
    unrounded = coeffs[np.newaxis, :] * scaling_factors[:, np.newaxis]
    absunrounded = np.abs(unrounded)
    rounded = np.round(unrounded)
    difference = np.abs(rounded - unrounded)
    difference[difference <= eps10] = 0.0
    num_under_threshold = np.sum(absunrounded < under_thresh_limit, axis=1)
    residuals = scaling_factors * np.nansum(difference / absunrounded, axis=1)

    max_threshold_limit = len(coeffs) - 1
    max_allowed_under_threshold = 0

    while max_allowed_under_threshold <= max_threshold_limit:
        # Filter based on the allowed threshold criteria
        valid_indices = num_under_threshold <= max_allowed_under_threshold

        if np.any(valid_indices):
            valid_scaling_factors = scaling_factors[valid_indices]
            valid_difference = difference[valid_indices]
            valid_unrounded = unrounded[valid_indices]

            # Calculate residuals for valid scaling factors
            residuals = valid_scaling_factors * np.nansum(valid_difference / np.abs(valid_unrounded), axis=1)

            # Find the scaling factor with the minimum residual
            min_residual_index = np.argmin(residuals)
            best_scaling = valid_scaling_factors[min_residual_index]

            return best_scaling

        # Increment max_allowed_under_threshold if no valid indices are found
        max_allowed_under_threshold += 1

    # If no valid scaling factor is found within the allowed threshold limit
    return None

##################################################################################################################
# Smoothing with 0th..5th order polynomials
##################################################################################################################

# Generate Savitzky-Golay coefficients and scaled results for various window sizes
# window_sizes = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63]
window_sizes = [3, 5]

##################################################################################################################

polyorder = 1
deriv = 0
delta = 1.0
use = 'conv'
linearSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 1
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    linearSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _linearSmooth[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(linearSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 2
deriv = 0
delta = 1.0
use = 'conv'
quadraticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 1
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    quadraticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quadraticSmooth[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(quadraticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 3
deriv = 0
delta = 1.0
use = 'conv'
cubicSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 1
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    cubicSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _cubicSmooth[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(cubicSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 4
deriv = 0
delta = 1.0
use = 'conv'
quarticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 1
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    quarticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quarticSmooth[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(quarticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 5
deriv = 0
delta = 1.0
use = 'conv'
quinticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 1
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    quinticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quinticSmooth[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(quinticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################
# 1st Drivative
##################################################################################################################

polyorder = 1
deriv = 1
delta = 1.0
use = 'conv'
linearSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    linearSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _linearDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(linearSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 2
deriv = 1
delta = 1.0
use = 'conv'
singQuadraticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    singQuadraticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quadraticFirstDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(singQuadraticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 3
deriv = 1
delta = 1.0
use = 'conv'
quadCubicSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    quadCubicSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _cubicFirstDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(quadCubicSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 4
deriv = 1
delta = 1.0
use = 'conv'
cubeQuarticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    cubeQuarticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quarticFirstDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(cubeQuarticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 5
deriv = 1
delta = 1.0
use = 'conv'
quartQuinticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    quartQuinticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quinticFirstDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(quartQuinticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################
# 2nd Drivative
##################################################################################################################

polyorder = 2
deriv = 2
delta = 1.0
use = 'conv'
singQuadraticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    singQuadraticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quadraticSecondDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(singQuadraticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 3
deriv = 2
delta = 1.0
use = 'conv'
quadCubicSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    quadCubicSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _cubicSecondDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(quadCubicSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 4
deriv = 2
delta = 1.0
use = 'conv'
cubeQuarticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    cubeQuarticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quarticSecondDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(cubeQuarticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 5
deriv = 2
delta = 1.0
use = 'conv'
quartQuinticSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = findScale(reduced_coeffs, max_s)
    else:
        coeffs = np.zeros(window_length)
        coeffs[pos] = 0
        reduced_coeffs = coeffs[:window_length // 2 + 1]
        scale = 1
    scaled_coeffs = np.round(reduced_coeffs * scale).astype(int).tolist()
    # Ensure the length of the coefficients list is max_coeffs_len, padding with zeros if necessary
    scaled_coeffs += [0] * (max_coeffs_len - len(scaled_coeffs))
    scaled_coeffs.insert(0, scale)
    quartQuinticSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"int32_t _quinticSecondDerivative[{len(window_sizes)}][{max_coeffs_len+1}] = {{")
for i, coeffs in enumerate(quartQuinticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

