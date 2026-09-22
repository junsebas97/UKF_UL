# Unscented Kalman Filter with Uncertain Inputs (UKF-UL)

This repository implements the UKF-UL [1]. This filter estimates states and uncertain forces using Gaussian random walks, the unscented transform, and single-step time-integration algorithms.

The repository contains a [numerical validation and comparison with the AUKF and the UKF-UI [2]](./numerical_validation.m), a [sensitivity analysis on the noise level](./sensitivity_measurement_noise.m), and a [sensitivity analysis on the number of measurements](./sensitivity_number_measurement.m). In the examples, the analysed system is a 6-DOFs shear chain with hysteretic springs governed by the Bouc-Wen model.

## Bibliography:
[1] J. S. Delgado Trujillo, J. Tott-Buswell, S. Jalbi, J. Hilton, M. Pandey, L. J. Prendergast, A nonlinear Bayesian filter for structural systems subjected to uncertain loads, Journal of Sound and Vibration. [In press]

[2] Y. Lei, D. Xia, K. Erazo, S. Nagarajaiah, A novel unscented Kalman filter for recursive state-input-system identification of nonlinear systems, Mechanical Systems and Signal Processing 127 (2019) 120–135. doi:10.1016/j.ymssp.2019.03.013.