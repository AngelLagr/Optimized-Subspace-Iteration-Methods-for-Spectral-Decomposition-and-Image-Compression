# Optimized Subspace Iteration Methods for Spectral Decomposition and Image Compression

This repository contains the implementation of advanced **Subspace Iteration Methods** for spectral decomposition of matrices, with a focus on applying these methods to **image compression**. 

## Project Overview
We explore various versions of the **Subspace Iteration Method**, starting from the **Power Method** and extending to optimized versions including **Rayleigh Projection** and acceleration techniques. The project includes a performance comparison of custom algorithms against MATLAB’s built-in `eig` function for eigenvalue computations, particularly applied to image compression tasks.

## Key Features
- **Subspace Iteration**: Various implementations for improved spectral decomposition.
- **Rayleigh Projection**: Faster convergence for eigenvalue computation.
- **Image Compression**: Using matrix decomposition for efficient image compression.
- **Performance Metrics**: Comparative analysis with MATLAB’s `eig` function.

## Methods Implemented
1. **Power Method**: A basic iterative method for eigenvalue computation.
2. **Subspace Iteration (v1)**: A refined method without projection.
3. **Subspace Iteration with Rayleigh Projection (v2)**: Enhanced with Rayleigh projection.
4. **Subspace Iteration with Rayleigh Projection and Acceleration (v3)**: Further optimized for large matrices.

## Results
- The **MATLAB `eig` function** is generally the fastest and most efficient.
- The **accelerated subspace iteration** method performs well for large matrices in image compression applications.
