# Analytical and deep-learning assisted reconstruction code for 3I-SIM

This is the supplementary software and source code for our manuscript "Triangle beam interference structured illumination microscopy", which is accepted by Nature Photonics 🚩.

This work introduces a new type of structured illumination microscopy that uses radially polarized three-beam interference.


Currently, the package contains the .m source code and a .exe files containing GUI of the 3I-SIM reconstruction code, and .py files for deep learning 3I-SIM reconstruction.

# System requirements
The software is implemented with MATALB 2020b and tested in Windows 10&11 environments.


# Installation guide
The package does not need additional installation steps. 

To run the conventional reconstruction algorithm with the .m file, it is only required to install MATLAB at any version (this code is written based on MATLAB 2020b, different MATLAB version may have different function rules that influence the functioning).

To run the .exe files, the MATLAB software is not necessary but one should download and install the 15 MATLAB runtime version 2020b for compiling from: https://ww2.mathworks.cn/products/compiler/matlab-runtime.html.

# Demo
Some test data and parameters are attached with the project for demonstration. 

# Running time
It usually takes 5 to 20 seconds to reconstruct a 1024x1024 super-resoltuion image from a 512x512x7 image stack, depending on the CPU performance. 

# Instructions for use
A detailed UserManual is attached with the project. For any questions in usage, please contact: houyiwei@stu.pku.edu.cn
