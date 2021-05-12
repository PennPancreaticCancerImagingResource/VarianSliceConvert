# VarianSliceConvert
Convert fdf file to nifti (or other) format

## Install
CMake and ITK are needed to compile this program

Create a build directory, cd into it, and run: ccmake YOUR_PATH/VarianSliceConvert

If cmake doesn't find ITK, you will need to point it to your install location. Then hit 'c' for configure, then 'g' for generate. Finally run 'make' and the executable should be compiled


## Usage
VarianSliceConvert input.fdf output.nii

In addition to nifti (shown above), any format supported by the default IO may be used for output.
