# isofin
Isogeometric Fiber Networks Simulation tool. A complete pipeline to pre-process, solve and post-process random fiber networks using Isogeometric Analysis.

![](logo/ChatGPT_logo.webp)

## Requirements

* Matlab: pre-processing and data analysis
* Paraview: post-processing and visualization
* C++ Standard Libraries
* OpenMP: library for parallel programming, required by our solver

## Overview

![](logo/Pipeline.png)

The following is a guide to using Isofin. MATLAB scripts are used to generate networks, input files, and mesh files. The input files and mesh files are used in the C++ code to run a nonlinear finite element simulation with isogeometric beams. Then, the results output from this C++ code are analyzed using MATLAB scripts, and converted into a form readable by Paraview. Detailed steps are given below, along with relevant information.

## Creating a fiber network
1.	Navigate to “isofin/Matlab/Network Generation”. “Fractal_NetworkGenerator_3D.m” will create a 3D fiber network of straight fibers. Further steps are required to add undulations to the network. Open “Fractal_NetworkGenerator_3D.m” and set the network parameters to generate your desired network. 
a.	Here, you have control over the branching angles, fiber length, number of networks to create, and the number of fibers your networks should have.
b.	Note: we advise against setting NumIter too high, as it will run for a significant amount of time due to the recursive nature of “Branch_It.m”. The best value for this will depend on your machine. We recommend starting with NumIter = 10 and adjusting from there.
2.	If you are using a straight fiber network, proceed to the next step. If you wish to use an undulated fiber network, you will need to use “Undulating_Networks.m”. Here, you may adjust the period to modify the undulations of the network. Be sure to properly set num_fibers and Case to match the network you want to add undulations to.

## Creating input and mesh files
1.	Navigate to “isofin/Matlab/Input File Generation”. The file “Input_File_Generation.m” will turn the networks created in the previous step into input and mesh files suitable for analysis in our C++ code. Here, you will set the type of network you are using, deformation mode of the network, and which fiber network you will use. 
a.	In regards to the network type, our scripts can handle 2D and 3D networks, but you will need to generate your own Voronoi networks to generate 2D input and mesh files.
b.	Within this script, you have control over the element size, the order of NURBS functions, and network parameters for the trimming of fibers and insertion of control points for the bending strips.
c.	We have built into the code uniaxial stretching, biaxial stretching, and simple shear. However, you may adjust the BCs as desired. You can apply displacement or force BCs. You can control the number of threads on which to run the simulation, the number of increments you want, the max number of iterations per increment, the max number of attempts per increment, and convergence criterion.
d.	Make note of where you save your input and mesh files. We opt for saving them in “isofin/Input_files” and “isofin/Mesh_files”. You will need to adjust the paths in the C++ codes to find the files if you want to organize your files differently.

## Running simulations  
1.	Within “isofin/Main” are three C++ files. “RN_2D_Trial_UATX.cpp” is for any 2D networks you generate, while “RN_3D_Debug.cpp” and “RN_3D_Trial_UATX.cpp” are for any 3D networks you generate. To compile the code, please use the following command using “RN_3D_Trial_UATX.cpp” as an example:
a.	 “<compiler> RN_3D_Trial_UATX.cpp -o Analysis_3D -O3 -fopenmp -std=c++17”
b.	We recommend using a C++ standard no newer than C++17, as newer standards will cause the Eigen library to generate warnings. We cannot guarantee convergence of the simulations with newer C++ standards either. 
c.	The results files will be stored in “isofin/Results/Outputs/3D” (or 2D if running a 2D network), unless otherwise changed. This is the directory in which the scripts in “isofin/Matlab/Output Generation” will look for the results.

## Generating Paraview compatible output files
1.	In “isofin/Matlab/Output Generation”, there are two scripts, one for generating *vtk files for straight networks, and one for undulated networks. They are named accordingly, “Paraview_VTK_Generator_Line_3D_Straight.m” generated *vtk’s for straight fiber networks, and “Paraview_VTK_Generator_Line_3D_Undulated.m” generate *vtk’s for undulated fiber networks. Note that a *vtk file will be generated for each increment in the simulation. The script also calculates the energy and energy ratios for the fibers in the network. After running one of these scripts, you will have *vtk files for your simulation, and may perform further analyses on your results.

Thank you for using Isofin! If you have any questions or comments regarding this project, please reach out to us via email.
The original author of these codes is Soham Mane (sohammane@utexas.edu). Both Matthew J Lohr (mlohr@utexas.edu) and Sotiris Kakaletsis (kakalets@utexas.edu) edited the codes, prepared them for sharing, and created this documentation.
