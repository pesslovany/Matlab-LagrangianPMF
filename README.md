# Matlab-LagrangianPMF

This repository implements Lagrangian grid based filter, Particle filter and Unscented Kalman Filter. The filters are used for position estimation of a vehicle based on a real world barometric measurements. The code is meant to complement a publication: TBD DOI. For more information please read the paper. 

__If you are using this code for a publication, please cite the paper.__

The parameters of the filters and models used (where UKF can be switched off) are in the file initModel.m. 

The code is run from the main.m. Matlab data data_UKFdivergence_3D and data_UKFdivergence_4D can be loaded from the main.m file lines 39/40, they showcase data for which UKF fails. There are two mex files, they should be compiled and ready to use with windows, intel mac and silicone mac. If there is an issue with these files, please comment out the parts that call them and change the lines for the commented out parts that does not need them. That is

For binarySearch: 
```
% Use this in case the mex file not working
% for ind = 1:1:noParts
%     I(ind) = find(cumW >= randomN(ind),1, 'first');
% end
I = binarySearch(cumW, randomN, 'first');
``` 
For inplaceprod:
```
% Use this in case the mex file not working
% mesPdfDotDeltasTens = mesPdfDotDeltasTens.*convKerTens;
% Convolution in frequency domain (inplaceprod saves memory)
inplaceprod(mesPdfDotDeltasTens, convKerTens);
```

## User can swich between number of models:

### 4D Model with 1D/3D measurement

This model can be used by modelChoose = 4;/modelChoose = 41; (3D measurement 1D measurement)
 
The state is the vehicle horizontal position $\mathbf{p}_k^\mathrm{W} = [p_k^{x,\mathrm{W}}, p_k^{y,\mathrm{W}}]^T$ [m] and velocity $\mathbf{v}_k^\mathrm{W} = [v_k^{x,\mathrm{W}}, v_k^{y,\mathrm{W}}]^T$ [m/sec] in a world (W) frame aligned with the geographic north, while the measurement is given by the barometric altimeter, which reads the vehicle altitude above the mean sea level (MSL) $\hslash_k^\mathrm{MSL}$ [m], and the odometer, which provides the vehicle velocity in the *body* (B) frame $\mathbf{v}_k^\mathrm{B} = [v_k^{x,\mathrm{B}}, v_k^{y,\mathrm{B}}]^T$ [m/sec].

We assume that the heading of the vehicle is aligned with $\mathbf{v}_k^\mathrm{B}$; consequently, the W and B frames are rotated by the heading angle $\psi_k$ (angle between $\mathbf{v}_k^\mathrm{B}$ and $\mathbf{v}_k^\mathrm{W}$) and the respective direction cosine matrix (DCM) is 

$$
\mathbf{C}_k = \begin{bmatrix}
    \cos(\psi_k) & -\sin(\psi_k) \\
    \sin(\psi_k) & \cos(\psi_k)
\end{bmatrix}.
$$

The model reads:

$$
\mathbf{x}_{k+1} = \begin{bmatrix}
    \mathbf{p}_k^{\mathrm{W}} \\
    \mathbf{v}_k^\mathrm{W} 
\end{bmatrix} = \begin{bmatrix}
    1 & 0 & T & 0 \\
    0 & 1 & 0 & T \\
    0 & 0 & 1 & 0 \\
    0 & 0 & 0 & 1
\end{bmatrix} \mathbf{x}_k + \mathbf{w}_k,
$$

$$
\mathbf{z}_k = \begin{bmatrix}
    \hslash_k^\mathrm{MSL} \\
    \mathbf{v}_k^\mathrm{B} 
\end{bmatrix} = \begin{bmatrix}
    \mathrm{terMap}(\mathbf{p}_k^\mathrm{W}) \\
    \mathbf{C}_k \mathbf{v}_k^\mathrm{W} 
\end{bmatrix} + \mathbf{v}_k.
$$

For 1D variant the measurement is only the vehicle altitude above mean sea level. The 1D variant does not offer enough information about the vehicle position and all filters diverge when the vehicle gets to a place where the map is not informative enought i.e. the terrain is flat. 

### 3D Model

The state is the vehicle horizontal position $\mathbf{p}_k^\mathrm{W} = [p_k^{x,\mathrm{W}}, p_k^{y,\mathrm{W}}]^T$ [m] and barometer bias $\mathbf{b}_k$ [m]. The measurement is given by the barometric altimeter, which reads the vehicle altitude above the mean sea level (MSL) $\hslash_k^\mathrm{MSL}$ [m]. Input $\mathbf{u}_k$ is given by shift, this shift can be given by a combination of odometer (information in the body frame) and magnetometer (to tie the odometer information to the world frame). For our purpose the odometer/magnetometer joined information is simulated.

The model reads:

$$
\mathbf{x}_{k+1} = \begin{bmatrix} \mathbf{p}_k^{\mathrm{W}} \\ \mathbf{b}_k
\end{bmatrix} = \mathbf{x}_k + \mathbf{u}_k + \mathbf{w}_k,
$$

$$
\mathbf{z}_k = \hslash_k^\mathrm{MSL}= \mathrm{terMap}(\mathbf{p}_k^\mathrm{W}) + \mathbf{v}_k.
$$
