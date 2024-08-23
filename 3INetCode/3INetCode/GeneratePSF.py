import numpy as np
from math import exp
import torch

def GPSF(lambda_em, pixelsize, NA):
    w = 9
    # Create mesh of X and Y coordinates
    X, Y = np.meshgrid(np.linspace(0, w - 1, w), np.linspace(0, w - 1, w))
    X = X - w / 2
    Y = Y - w / 2
    R = X ** 2 + Y ** 2
    # Calculate sigmap using the Bessel functions and other terms
    lambda_em /= pixelsize
    sigmap = 0.22 * lambda_em / NA

    # Calculate the PSF using the exponential function and fftshift for centering
    PSF = np.exp(-R / (
                2 * sigmap ** 2))  # Normalization is done here instead of at the end due to Python's behavior with division and exponentiation order of operations
    PSF = PSF / np.sum(np.sum(PSF))
    PSF = PSF.reshape(1, 1, 9, 9)
    PSF = torch.tensor(PSF, dtype=torch.float32)
    return PSF