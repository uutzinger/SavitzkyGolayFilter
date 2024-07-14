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

# We have 32bit intger and floating point math on ESP32 and ESP8266
# We have 64bit intger and floating point math on teensy
# 
# Scaling and filter coefficients should be limitted to minimize nummerical overflow.
# If data is 16 bit and we convolve window_size number of data points we can set a conservative limit for scaling factor as shown below:

max_s = np.iinfo(np.int16).max * 25
# max_s = 10240000 # Sacitzky-Golay have in their paper scaling coeffcients up to 430 Million,

def findScale(coeffs, max_scale=256000):
    # Searching for best scaling factor is vectorized and  will calcualte all residuals for all scaling factors at once
    # Simplified explanation:
    # 1) error = round(coeffs*scaling) / scaling - coeffs
    # 2) residual = sum(abs(error))
    # 3) find smallest residual which is best scaling factor
    eps = np.finfo(np.float64).eps # machine precision
    scaling_factors = np.arange(1, max_scale) # all scaling factors
    rounded_coeffs_scaled = np.round(coeffs[np.newaxis, :] * scaling_factors[:, np.newaxis])
    errors = np.abs(rounded_coeffs_scaled/scaling_factors[:, np.newaxis] - coeffs[np.newaxis, :])
    errors[errors <= eps] = 0.
    residuals = np.sum(errors, axis=1)
    min_residual_index = np.argmin(residuals)
    best_scaling = scaling_factors[min_residual_index]
    return best_scaling

######################################################################################################################
# SMALL MODE, window size up to 25
######################################################################################################################

window_sizes = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
eps = np.finfo(np.float64).eps

##################################################################################################################
# Smoothing with 1st..5th order polynomials
##################################################################################################################

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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
#print(f"std::vector<std::vector<int32_t>> _linearSmooth = {{")
print(f"std::vector<std::vector<int32_t>> _linearSmooth = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quadraticSmooth = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _cubicSmooth = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quarticSmooth = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quinticSmooth = {{")
for i, coeffs in enumerate(quinticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 6
deriv = 0
delta = 1.0
use = 'conv'
sexicSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        coeffs[np.abs(coeffs) <= eps] = 0.
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
    sexicSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"std::vector<std::vector<int32_t>> _sexicSmooth = {{")
for i, coeffs in enumerate(sexicSmooth):
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _linearFirstDerivative = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.

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
print(f"std::vector<std::vector<int32_t>> _quadraticFirstDerivative = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _cubicFirstDerivative = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quarticFirstDerivative = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quinticFirstDerivative = {{")
for i, coeffs in enumerate(quartQuinticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 6
deriv = 1
delta = 1.0
use = 'conv'
sexicSmooth = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        coeffs[np.abs(coeffs) <= eps] = 0.
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
    sexicSmooth.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"std::vector<std::vector<int32_t>> _sexicFirstDerivative = {{")
for i, coeffs in enumerate(sexicSmooth):
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quadraticSecondDerivative = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _cubicSecondDerivative = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quarticSecondDerivative = {{")
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
        coeffs[np.abs(coeffs) <= eps] = 0.
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
print(f"std::vector<std::vector<int32_t>> _quinticSecondDerivative = {{")
for i, coeffs in enumerate(quartQuinticSmooth):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

##################################################################################################################

polyorder = 6
deriv = 2
delta = 1.0
use = 'conv'
sexic = []

max_window_size = max(window_sizes)
max_coeffs_len = max_window_size // 2 + 1

for window_length in window_sizes:
    pos = window_length // 2
    if polyorder < window_length:
        coeffs = savgol_coeffs(window_length, polyorder, deriv=deriv, delta=delta, pos=pos, use=use)
        coeffs[np.abs(coeffs) <= eps] = 0.
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
    sexic.append(scaled_coeffs)

# Print the C-style array definition
print("")
print(f"std::vector<std::vector<int32_t>> _sexicSecondDerivative = {{")
for i, coeffs in enumerate(sexic):
    formatted_coeffs = ', '.join(f"{x: >4}" for x in coeffs)
    print(f"  {{ {formatted_coeffs} }}, // Window size {window_sizes[i]}")

print("};")

