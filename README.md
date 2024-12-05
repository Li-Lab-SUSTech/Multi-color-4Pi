# Multi-color-4Pi
 Ratiometric-4Pi is a graphics processing unit (GPU) based global fitting algorithm for 4Pi-SMLM with flexible PSF modeling and parameter sharing, to extract maximum information from 4Pi single molecule data and achieved both good color separation and optimal 3D resolution. By partially linking the photon parameters between channels with interference difference of π during global fitting of the multi-channel 4Pi single molecule data, we showed on simulated data that the loss of the localization precision is minimal compared with the theoretical minimum uncertainty, the Cramer-Rao lower bound (CRLB). Our algorithm is implemented in GPU and the fitting speeds is more than 38 times faster than the CPU based code.
 
 <img src="https://github.com/Li-Lab-SUSTech/Ratiometric-4Pi/blob/main/Figure/Fig_1_Schematic of 4Pi-SMLM_1.png" width = 80%  alt="workflow overview" align=center />

This code comes with the paper: "[Ratiometric 4Pi single-molecule localization with optimal resolution and color assignment](https://opg.optica.org/ol/abstract.cfm?uri=ol-47-2-325)".

If you use this code for your research, please cite our paper:

* Jianwei Chen, Haoyu Wang, Zhaojun Lin, and Yiming Li, "Multicolor 4Pi single molecule localization based on differences in interference patterns," Opt. Lett. 47, 325-328 (2022)

# Requirements
Matlab R2019a or newer  

The GPU fitter requires:
  
  - Microsoft Windows 7 or newer, 64-bit
  - CUDA capable graphics card, minimum Compute Capability 6.1
  - CUDA 11.3 compatible graphics driver (for GeForce products 471.41 or later)

The CPU version runs on macOS and Microsoft Windows 7 or newer, 64-bit

 # How to run
 Examples code are avalible in file **dual_color_separation_for_random_NP.m, dual_color_separation_for_uniform_NP.m, three_color_separation_for_random_NP.m, three_color_separation_for_uniform_NP.m, Compare_CRLB_with_Saved_and_Ratio_4Pi.m**. 
 
# Contact
For any questions / comments about this software, please contact [Li Lab](https://faculty.sustech.edu.cn/liym2019/en/).

# Copyright 
Copyright (c) 2024 Li Lab, Southern University of Science and Technology, Shenzhen.
