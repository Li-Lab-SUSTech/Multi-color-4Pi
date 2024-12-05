# Multi-color-4Pi
 Multi-color-4Pi is a GPU-based multicolor algorithm for 4Pi-SMLM that extracts color information directly from 4Pi single-molecule data by utilizing wavelength-dependent point spread functions (PSFs). By globally fitting data from multiple interference channels, this method achieves high-precision localization and excellent color separation for single molecules of different colors. The algorithm requires no hardware modifications and can be seamlessly integrated with existing 4Pi systems.
 
 <img src="https://github.com/Li-Lab-SUSTech/Multi-color-4Pi/blob/main/Figure/Fig1_v2.png" width = 60%  alt="workflow overview" align=center />

<!This code comes with the paper: "[Ratiometric 4Pi single-molecule localization with optimal resolution and color assignment](https://opg.optica.org/ol/abstract.cfm?uri=ol-47-2-325)".

If you use this code for your research, please cite our paper:

* Jianwei Chen, Haoyu Wang, Zhaojun Lin, and Yiming Li, "Multicolor 4Pi single molecule localization based on differences in interference patterns," Opt. Lett. 47, 325-328 (2022)>

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
