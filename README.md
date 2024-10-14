# Matlab-LagrangianPMF

This repository implements Lagrangian grid based filter for Terrain aided navigation, and compares it with other state of the art filters. The filters are used for position estimation of a vehicle based on a real world barometric measurements. The code is meant to complement a publication: TBD DOI. For more information please read the paper. 

__If you are using this code for a publication, please cite the paper.__

The parameters of the filters and models used for the publication results are in the file initModel.m. 

The Lagrangiang grid based filter code is run from the main.m. The comparison with other state of the art filters is in the folder comparison, and can be run by mainComparison.m. There are two mex files, they should be compiled and ready to use with windows, linux, intel mac and silicone mac. If there is an issue with these files, please comment out the parts that call them and change the lines for the commented out parts that does not need them. That is

For binarySearch in mainComparison: 
```
% Use this in case the mex file not working
% for ind = 1:1:noParts
%     I(ind) = find(cumW >= randomN(ind),1, 'first');
% end
I = binarySearch(cumW, randomN, 'first');
``` 
For inplaceprod usade by Lagrangian grid based filter:
```
% Use this in case the mex file not working
% mesPdfDotDeltasTens = mesPdfDotDeltasTens.*convKerTens;
% Convolution in frequency domain (inplaceprod saves memory)
inplaceprod(mesPdfDotDeltasTens, convKerTens);
```

## User can swich between number of models:

### 4D Model
 
# Model Description

The state consists of the vehicle's horizontal position $\mathbf{p}_k^\mathrm{W} = [p_k^{x,\mathrm{W}}, p_k^{y,\mathrm{W}}]^T$ [m] and velocity $\mathbf{v}_k^\mathrm{W} = [v_k^{x,\mathrm{W}}, v_k^{y,\mathrm{W}}]^T$ [m/s] in a world frame (W) aligned with the geographic north. The total state dimension is $n_x = 4$, consisting of two components each for position and velocity. The measurement vector, with dimension $n_z = 3$, consists of the barometric altimeter reading $\hslash_k^\mathrm{MSL}$ [m] and the velocity measured by the odometer in the body (B) frame $\mathbf{v}_k^\mathrm{B} = [v_k^{x,\mathrm{B}}, v_k^{y,\mathrm{B}}]^T$ [m/s].

The state transition is modeled with a time step $T = 1$ [s] and includes process noise characterized by a covariance matrix $\mathbf{Q}$. The discrete-time system dynamics are described by the following equation:

$$
\mathbf{x}_{k+1} = \mathbf{F} \mathbf{x}_k + \mathbf{w}_k,
$$
where the system matrix $\mathbf{F}$ is defined as:

$$
\mathbf{F} = \begin{bmatrix}
    1 & 0 & T & 0 \\
    0 & 1 & 0 & T \\
    0 & 0 & 1 & 0 \\
    0 & 0 & 0 & 1
\end{bmatrix},
$$
and $\mathbf{w}_k$ is the process noise with covariance matrix $\mathbf{Q}$ given by:

$$
\mathbf{Q} = q \begin{bmatrix}
    \frac{T^3}{3} & 0 & \frac{T^2}{2} & 0 \\
    0 & \frac{T^3}{3} & 0 & \frac{T^2}{2} \\
    \frac{T^2}{2} & 0 & T & 0 \\
    0 & \frac{T^2}{2} & 0 & T
\end{bmatrix},
$$
where $q = 1$ [m$^2$/s$^3$] is a noise parameter.

The measurements are modeled as:

$$
\mathbf{z}_k = \begin{bmatrix}
    \hslash_k^\mathrm{MSL} \\
    \mathbf{v}_k^\mathrm{B}
\end{bmatrix} = \begin{bmatrix}
    \text{terMap}(\mathbf{p}_k^\mathrm{W}) \\
    \mathbf{C}_k \mathbf{v}_k^\mathrm{W}
\end{bmatrix} + \mathbf{v}_k,
$$
where $\mathbf{v}_k$ is the measurement noise with covariance matrix:

$$
\mathbf{R} = \text{diag}([3, 1, 1]^2),
$$
representing the noise variances of the altimeter and odometer measurements.

## Particle Filter Parameters

- The model employs $N_\mathrm{pa} = [51, 51, 41, 41]$ grid points for each state axis, leading to a total of $N = \prod N_\mathrm{pa}$ points. The number of particles used in the Particle Filter is $\text{noPart} = 1.8N$.
- The scaling factor for the grid, which determines how many standard deviations the grid spans, is set to $s_\mathrm{factor} = 5$.
- The effective sample size threshold is $\mathrm{essThrd} = \frac{2}{3} \text{noPart}$.

## Initial Conditions

- The initial state mean $\mathbf{x}_0$ is set according to GNSS coordinates: $\mathbf{x}_0 = [p_0^{x,\mathrm{W}}, p_0^{y,\mathrm{W}}, v_0^{x,\mathrm{W}}, v_0^{y,\mathrm{W}}]^T$.
- The initial state covariance matrix $\mathbf{P}_0$ is:

$$
\mathbf{P}_0 = \begin{bmatrix}
    25 & 0 & 0 & 0 \\
    0 & 25 & 0 & 0 \\
    0 & 0 & 0.5 & 0 \\
    0 & 0 & 0 & 0.5
\end{bmatrix}.
$$

## Measurement Model

The measurement equation is defined as:

$$
\mathbf{h}(\mathbf{x}_k, \mathbf{v}_k, k) = \begin{bmatrix}
    \text{terMap}(\mathbf{p}_k^\mathrm{W}) \\
    \mathbf{C}_k \mathbf{v}_k^\mathrm{W}
\end{bmatrix} + \mathbf{v}_k,
$$
where the measurement noise for the odometer is added after generating the measurement using the measurement noise covariance matrix $\mathbf{R}$.

### 3D Model

The state is the vehicle horizontal position $\mathbf{p}_k^\mathrm{W} = [p_k^{x,\mathrm{W}}, p_k^{y,\mathrm{W}}]^T$ [m] and barometer bias $\mathbf{b}_k$ [m]. The measurement is given by the barometric altimeter, which reads the vehicle altitude above the mean sea level (MSL) $\hslash_k^\mathrm{MSL}$ [m]. Input $\mathbf{u}_k$ is given by shift; this shift can be given by a combination of odometer (information in the body frame) and magnetometer (to tie the odometer information to the world frame). For our purpose, the odometer/magnetometer joined information is simulated.

\subsection*{System Parameters}
The model parameters are defined as follows:
\begin{itemize}
    \item State dimension: $n_x = 3$
    \item Measurement dimension: $n_z = 1$
    \item Time step: $\Delta t = 1$ [s]
    \item System noise covariance: 
    $$ \mathbf{Q} = \mathrm{blkdiag}(36 \cdot \mathbf{I}_2, 0.05) $$
    \item Input: 
    $$ \mathbf{u}_k = \begin{bmatrix} [\mathrm{diff}(\text{souradniceGNSS}(1:2, \text{timeSteps}), 1, 2) \; [0; 0]] + \sqrt{\mathbf{Q}(1:2, 1:2)} \cdot \mathbf{randn}(2, \text{endTime}) \\ 0 \end{bmatrix} $$
\end{itemize}

\subsection*{PMF Parameters}
The PMF parameters are as follows:
\begin{itemize}
    \item Number of points per axis: $\mathbf{Npa} = [71, 71, 21]$
    \item Total number of points: $N = \prod(\mathbf{Npa})$
    \item Number of particles for multinomial resampling PF: $N_{part} = 250000$
    \item Effective sample size threshold for PF: $ess_{thrd} = \frac{2}{3} \cdot N_{part}$
    \item Scaling factor (number of sigmas covered by the grid): $s_{factor} = 5$
    \item Mean values of components of measurement noise: $\mathbf{meanV} = 0$
    \item Weights: $w_V = 1$
\end{itemize}

\subsection*{Initial Conditions}
The initial conditions are defined as:
\begin{itemize}
    \item Initial state: 
    $$ \mathbf{x} = \begin{bmatrix} \text{souradniceGNSS}(1:2, \text{timeSteps}) \\ h_B - \text{vysky}(\text{souradniceGNSS}(1, \text{timeSteps}), \text{souradniceGNSS}(2, \text{timeSteps})) \end{bmatrix} $$
    \item Initial condition mean: 
    $$ \mathbf{meanX0} = \begin{bmatrix} \text{souradniceGNSS}(1:2, \text{timeSteps}(1)) \\ 0 \end{bmatrix} $$
    \item Initial condition variance: 
    $$ \mathbf{varX0} = \begin{bmatrix} 25 & 0 & 0 \\ 0 & 25 & 0 \\ 0 & 0 & 10 \end{bmatrix} $$
\end{itemize}

\subsection*{Dynamics}
The model dynamics are represented as:
$$ \mathbf{F} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0.96925 \end{bmatrix} $$

\subsection*{Measurement Generation}
The measurement equation is given by:
$$ \mathbf{z} = \text{hfunct}(\mathbf{x}, 0, 0) + \sqrt{\mathbf{R}} \cdot \mathbf{randn}(n_z, \text{length}(\mathbf{x})) $$
where the measurements are updated as:
$$ \mathbf{z}(1, :) = h_B $$

\subsection*{UKF Parameters}
The UKF parameters are defined as:
\begin{itemize}
    \item $\kappa = 1$
    \item Number of sigma points: $SP_{num} = 2n_x + 1$
\end{itemize}

