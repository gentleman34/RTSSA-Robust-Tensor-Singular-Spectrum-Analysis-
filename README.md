# Robust Low-Rank Reconstruction of 3-D Seismic Data

This repository provides a MATLAB implementation of a robust low-rank
reconstruction method for 3-D seismic data with randomly missing samples.

The implementation combines discrete cosine transform (DCT), frequency-space
Hankel matrix construction, robust low-rank matrix factorization, Tukey
bisquare weighting, and alternating weighted least-squares optimization.

A synthetic 3-D seismic dataset is provided to demonstrate the reconstruction
procedure.

## Repository Structure

```text
RTF/
├── README.md
├── LICENSE
├── demo.m
├── data/
│   └── ori2.mat
└── functions/
    ├── drr_plot3d.m
    └── other_required_functions.m
```

The functions `customFunction` and `backSubstitution` are included as local
functions at the end of `demo.m`.

## Requirements

The code is implemented in MATLAB.

Required MATLAB functions include:

- `dct`
- `idct`
- `fft`
- `ifft`
- `svd`
- `qr`

Please make sure that all visualization functions used by `demo.m`, such as
`drr_plot3d`, are included in the MATLAB search path.

## Input Data

The synthetic seismic dataset is stored in:

```text
data/ori2.mat
```

The MAT file contains the variable:

```matlab
ori2
```

which represents the original complete 3-D seismic data.

In the current experiment, the data dimensions are:

```text
500 × 100 × 80
```

corresponding to a three-dimensional seismic volume.

The data are loaded using:

```matlab
load(fullfile('data', 'ori2.mat'));
data = ori2;
```

## Generation of Missing Data

Randomly missing samples are generated from the original seismic data using:

```matlab
missing_ratio = 0.1;
mask = rand(size(data)) > missing_ratio;
data_missing = data .* mask;
```

The default missing ratio is therefore:

```text
10%
```

The variable `data_missing` is the observed incomplete seismic tensor used as
the input to the reconstruction procedure.

## Method Overview

The reconstruction procedure consists of the following steps.

### 1. Discrete Cosine Transform

A DCT is first applied along the first dimension of the incomplete 3-D seismic
tensor:

```matlab
fZ = dct(data_missing, [], 1);
```

### 2. DCT-Slice Processing

The first 80 DCT-domain slices are processed independently:

```matlab
for hi = 1:80
```

For each slice, a two-dimensional matrix is extracted and transposed:

```matlab
DATA = squeeze(fZ(hi,:,:));
DATA = permute(DATA,[2,1]);
```

### 3. Frequency-Space Transformation

Each selected DCT slice is transformed into the frequency-space domain using
the fast Fourier transform (FFT).

The default parameters are:

```matlab
dt    = 0.005;
flow  = 0;
fhigh = 130;
```

Only frequency components between `flow` and `fhigh` are processed.

### 4. Hankel Matrix Construction

For each frequency component, a level-1 Hankel matrix is constructed from the
spatial samples.

The Hankel matrix dimensions are determined by:

```matlab
Lcol = floor(nx/2) + 1;
Lrow = nx - Lcol + 1;
```

### 5. Low-Rank Initialization

The Hankel matrix is initialized using truncated singular value decomposition
(SVD):

```matlab
[initial_u, initial_d, initial_v] = svd(datam);
```

The default rank is:

```matlab
k = 5;
```

The low-rank Hankel matrix is represented as

```text
M ≈ U V^T
```

where `U` and `V` are the two low-rank factors.

### 6. Robust Weight Estimation

The residual is calculated as:

```matlab
RL = datam - U * V';
```

The residual scale is estimated using the median absolute deviation (MAD):

```matlab
seigma = 1.4862 * median(abs(r - median(r)));
```

Tukey bisquare weights are then calculated from the normalized residuals.

The default Tukey parameter is:

```matlab
alpha = 5;
```

Large residuals receive small or zero weights, reducing their influence on the
low-rank reconstruction.

### 7. Alternating Weighted Least-Squares Optimization

The two low-rank factors `U` and `V` are updated alternately by solving weighted
least-squares problems.

QR factorization and back substitution are used to solve the corresponding
linear systems.

The current implementation uses:

```text
3 outer iterations
3 inner iterations
```

### 8. Hankel Matrix Reconstruction

After the alternating optimization, the reconstructed Hankel matrix is obtained
as:

```matlab
result = U * V';
```

Anti-diagonal averaging is then used to transform the reconstructed Hankel
matrix back to the corresponding frequency-space samples.

### 9. Inverse FFT

Conjugate symmetry is imposed on the frequency-domain data, followed by an
inverse FFT to reconstruct each processed DCT slice.

### 10. 3-D Reconstruction

The reconstructed DCT slices are stored in:

```matlab
fresult
```

Finally, the reconstructed 3-D seismic volume is obtained using the inverse DCT:

```matlab
construct_result = idct(fresult, [], 1);
```

Therefore, the main reconstructed output of this implementation is:

```matlab
construct_result
```

## Main Parameters

The default parameters used in `demo.m` are:

| Parameter | Default value | Description |
|---|---:|---|
| `missing_ratio` | 0.1 | Ratio of randomly missing samples |
| `alpha` | 5 | Tukey bisquare tuning parameter |
| `k` | 5 | Rank of the low-rank factorization |
| `dt` | 0.005 | Sampling interval |
| `flow` | 0 | Lower frequency bound |
| `fhigh` | 130 | Upper frequency bound |
| DCT slices | 1–80 | DCT-domain slices processed |
| Outer iterations | 3 | Number of robust weight updates |
| Inner iterations | 3 | Number of alternating factor updates |

These parameters can be modified directly in `demo.m`.

## Running the Code

Clone or download the repository and open MATLAB in the repository root
directory.

Run:

```matlab
main
```

The program will:

1. load the complete synthetic seismic dataset;
2. randomly remove 10% of the samples;
3. transform the incomplete data into the DCT domain;
4. perform robust low-rank reconstruction;
5. apply the inverse DCT;
6. generate the reconstructed 3-D seismic volume.

The final reconstructed data are stored in:

```matlab
construct_result
```

## Output

The main variables generated by the program are:

| Variable | Description |
|---|---|
| `data` | Original complete 3-D seismic data |
| `mask` | Random sampling mask |
| `data_missing` | Seismic data with randomly missing samples |
| `fZ` | DCT-domain representation of the incomplete data |
| `fresult` | Reconstructed DCT-domain tensor |
| `construct_result` | Final reconstructed 3-D seismic data |
| `elapsedTime` | Computational time of the reconstruction |

The reconstruction error can be visualized using:

```matlab
drr_plot3d(data - construct_result, [500,100,80], 0.2);
```

## Reproducibility

The repository provides the synthetic input data and MATLAB implementation
required to reproduce the reconstruction workflow.

Because the missing-data mask is generated randomly, different runs may produce
slightly different missing patterns and reconstruction results.

For exactly reproducible experiments, a fixed random seed can be specified
before generating the sampling mask, for example:

```matlab
rng(2026);
```

## Notes

The current implementation processes the first 80 DCT-domain slices.

The visualization function `drr_plot3d` is used only for displaying 3-D seismic
data and is not part of the core reconstruction algorithm.

## License

Please refer to the `LICENSE` file for licensing information.

## Citation

If you use this code in your research, please cite the corresponding paper.

Citation information will be updated after publication.
