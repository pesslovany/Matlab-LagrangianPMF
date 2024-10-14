# Matlab-LagrangianPMF

This repository implements Lagrangian grid based filter for Terrain aided navigation, and compares it with other state of the art filters. The filters are used for position estimation of a vehicle based on a real world barometric measurements. The code is meant to complement a publication: TBD DOI. For more information please read the paper. 

__If you are using this code for a publication, please cite the paper.__

The parameters of the filters and models used for the publication results are described here and also in the file initModel.m. 

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
 
#### Model Description

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

where $q = 1$ [m $^2$ / s $^3$ ] is a noise parameter.

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
\mathbf{R} = \text{diag}([9, 1, 1]),
$$

representing the noise variances of the altimeter and odometer measurements.

#### Filter Parameters

- The proposed LGbF employs $N_\mathrm{pa} = [51, 51, 41, 41]$ grid points for each state axis. The Particle Filter uses $1.8$ times number of points of the LGbF.
- The scaling factor for the grid, which determines how many standard deviations the LGbF grid spans, is set to $s_\mathrm{factor} = 5$.
- The effective sample size threshold is $\mathrm{essThrd} = \frac{2}{3}$ number of particles.
- The jitter covariance is set-up to be 1% of the measurement estimate covariance.
- The particular details of the Teixiera et. al. Particle filter were not provided in the paper therfore they have been set up so the filter works well: 80% of particles were drawn from prior, and 20% from the uniform distribution, the uniform distribution was set up to cover $s_\mathrm{factor} = 5$ as in the LGbF case.

#### Initial Conditions

- The initial state mean $\mathbf{x}_0$ is drawn from the GNSS ground truth
- The initial state covariance matrix $\mathbf{P}_0$ is:

$$
\mathbf{P}_0 = \begin{bmatrix}
    25 & 0 & 0 & 0 \\
    0 & 25 & 0 & 0 \\
    0 & 0 & 0.5 & 0 \\
    0 & 0 & 0 & 0.5
\end{bmatrix}.
$$




### 3D Model

The state is the vehicle horizontal position $\mathbf{p}_k^\mathrm{W} = [p_k^{x,\mathrm{W}}, p_k^{y,\mathrm{W}}]^T$ [m] and barometer bias $\mathbf{b}_k$ [m]. The measurement is given by the barometric altimeter, which reads the vehicle altitude above the mean sea level (MSL) $\hslash_k^\mathrm{MSL}$ [m]. Input $\mathbf{u}_k$ is given by shift; this shift can be given by a combination of odometer (information in the body frame) and magnetometer (to tie the odometer information to the world frame). For our purpose, the odometer/magnetometer joined information is simulated.

The model reads:

$$
\mathbf{x}_{k+1} = \begin{bmatrix} \mathbf{p}_k^{\mathrm{W}} \\ \mathbf{b}_k
\end{bmatrix} = \mathbf{x}_k + \mathbf{u}_k + \mathbf{w}_k,
$$

where the system matrix $\mathbf{F}$ is defined as:

$$
\mathbf{F} = \begin{bmatrix}
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 0 & 0.96925 \\
\end{bmatrix},
$$

The measurements are modeled as:

$$
\mathbf{z}_k = \hslash_k^\mathrm{MSL}= \mathrm{terMap}(\mathbf{p}_k^\mathrm{W}) + \mathbf{v}_k.
$$

the $\mathbf{w}_k$ is the process noise with covariance matrix $\mathbf{Q}$ given by:

$$
\mathbf{Q} = q \begin{bmatrix}
   36 & 0 & 0\\
   0 & 36 & 0\\
   0 & 0 & 0.05\\
\end{bmatrix},
$$

the $\mathbf{v}_k$ is the measurement noise with covariance matrix:

$$
R = 3,
$$

representing the noise variances of the altimeter.

#### Filter Parameters

- The proposed LGbF employs $N_\mathrm{pa} = [71, 71, 21]$ grid points for each state axis. The Particle Filter uses $250000$ particles.
- The scaling factor for the grid, which determines how many standard deviations the LGbF grid spans, is set to $s_\mathrm{factor} = 5$.
- The effective sample size threshold is $\mathrm{essThrd} = \frac{2}{3}$ number of particles.
- The jitter covariance is set-up to be 1% of the measurement estimate covariance.
- The particular details of the Teixiera et. al. Particle filter were not provided in the paper therfore they have been set up so the filter works well: 80% of particles were drawn from prior, and 20% from the uniform distribution, the uniform distribution was set up to cover $s_\mathrm{factor} = 5$ as in the LGbF case.

#### Initial Conditions

- The initial state mean $\mathbf{x}_0$ is drawn from the GNSS ground truth and map.
- The initial state covariance matrix $\mathbf{P}_0$ is:

$$
\mathbf{P}_0 = \begin{bmatrix}
    25 & 0 & 0  \\
    0 & 25 & 0 \\
    0 & 0 & 10 \\
\end{bmatrix}.
$$

