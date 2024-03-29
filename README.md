UKF TRACTOGRAPHY
================

About
-----

We present a framework which uses an unscented Kalman filter for performing
tractography. At each point on the fiber the most consistent direction is found
as a mixture of previous estimates and of the local model.

It is very easy to expand the framework and to implement new fiber representations 
for it. Currently it is possible to tract fibers using two different 1-, 2-, or 3-tensor 
methods. Both methods use a mixture of Gaussian tensors. One limits the diffusion 
ellipsoids to a cylindrical shape (the second and third eigenvalue are assumed to be 
identical) and the other one uses a full tensor representation.

__Authors__:
Yogesh Rathi (yogesh@bwh.harvard.edu), Stefan Lienhard, Yinpeng Li, Martin
Styner, Ipek Oguz, Yundi Shi, Christian Baumgartner (c.f.baumgartner@gmail.com)
Ryan Eckbo, Rinat Mukhometzianov

For references and further documentation, please see the [Slicer module homepage](https://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/UKFTractography).

Installation
------------

### 1. From Source

Checkout from github:

    git clone https://github.com/rmukh/ukftractography_dp.git

There are 3 ways to build this project from source, as a stand alone
superbuild, against a Slicer 4 build, and as a Slicer 4 extension build (which
is more of a test than anything).

#### a) Standalone Superbuild

    cd <build-dir>
    cmake <path-to-source> -DCMAKE_BUILD_TYPE=Release|Debug|RelWithDebInfo
    make
    make test

##### On Windows,
    cd <build-dir>
    cmake <path-to-source> -DCMAKE_BUILD_TYPE=Release|Debug|RelWithDebInfo
    cmake --build . --config Release|Debug|RelWithDebInfo
    cmake --build build -t test

#### b) Build with Slicer4

    cd <build-dir>
    cmake -DSlicer_DIR=<path-to-Slicer4-Superbuild>/Slicer-build <path-to-source>
    make
    make test

The compilation might not work with the latest version of compilers. The recommended are GCC of versions 5,6,7,8. You can install them in addition to your default system compilers if needed and specify in the cmake arguments as `-DCMAKE_C_COMPILER=gcc-V -DCMAKE_CXX_COMPILER=g++-V`, where **V** is the version number from the list of recommended ones.

**Note** that some of the tests are broken!
#### c) Build via Slicer ExtensionIndex build

Create local extension index following [these instructions](https://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/Build_ExtensionsIndex), containing at least `UKFTractography.s4ext` and `SlicerDMRI.s4ext` (required runtime dependency).

Notes:

* To manually test the `UKF Tractography` Slicer module, start Slicer using the launcher named `SlicerWithSlicerDMRI` available in `/path/to/SlicerDMRI-build/inner-build` directory. This ensure that the SlicerDMRI modules are loaded and that the required MRML diffusion nodes are registered (i.e vtkMRMLFiberBundleNode).

* It may be helpful to [test the exension upload](https://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/Build_ExtensionsIndex#Extension_build.2C_test.2C_package_and_upload_using_.60ExperimentalUpload.60_target) using your API key.

### 2. As a Slicer 4 Extension

Navigate to the Slicer Extension Manager and download `UKF Tractography` to
install it as a Slicer 4 module.

Basic Usage
-----------

### 1. As Command Line Module

The executable is called 'UKFTractography'. It can be found in:
    
    <build-dir>/UKFTractography-build/UKFTractography/bin/

The path may containt a subfolder after bin/ if you are using a superbuild with the choosen build type (Release, Debug, etc.).

In order to see all options run.

    ./UKFTractography --help 

In the source directory of the project you will find a shell script called 'sample_run.sh'
It should give you an idea of what a function call could look like. 

Files dataset_Mask.nrrd and seeds_full_cc.nrrd in Input folder are mask and seed files of subject 100307
in hcp dataset, download the subject's preprocessed diffusion MRI data from https://db.humanconnectome.org/ 

### 2. As Slicer 4 module

Navigate to the Slicer Extension Manager and download `UKF Tractography` to
install it as a Slicer 4 module.  There will be 3 modules under
`Diffusion-->Tractography`: `UKF Tractography`, `vtk2mask`, and `vtkFilter`.

Notes
-----------

Several steps in the SuperBuild process download additional git repositories as CMake external projects via `https://` protocol.
You can find some running examples in the examples folder.
