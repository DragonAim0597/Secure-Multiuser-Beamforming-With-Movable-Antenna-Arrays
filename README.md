# Secure Multiuser Beamforming With Movable Antenna Arrays

The code for the paper **Z. Cheng, B. Zhao, C. Ouyang, and X. Zhang, “Secure multiuser beamforming with movable antenna arrays,” Jan. 2026.** 


Abstract: A movable antennas (MAs)-enabled secure multiuser transmission framework is developed to enhance physical-layer security. Novel expressions are derived to characterize the achievable sum secrecy rate based on the secure channel coding theorem. On this basis, a joint optimization algorithm for digital beamforming and MA placement is proposed to maximize the sum secrecy rate via fractional programming and block coordinate descent. In each iteration, every variable admits either a closed-form update or a low-complexity one-dimensional or bisection search, which yields an efficient implementation. Numerical results demonstrate the effectiveness of the proposed method and show that the MA-enabled design achieves higher secrecy rates than conventional fixed-position antenna arrays.


## Running the simulations

### Prerequisites

- [MATLAB](https://uk.mathworks.com/products/matlab.html)

### Launch

Run `Compare_SNR.m` to plot Fig. 3 in this paper.

Run `Compare_Antenna.m` to plot Fig. 4(a) in this paper.
