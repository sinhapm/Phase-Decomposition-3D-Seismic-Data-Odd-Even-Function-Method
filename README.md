# Phase-Decomposition-3D-Seismic-Data-Odd-Even-Function-Method
The Program performs the 3D seismic data Phase Decomposition based on Odd-Even function. It computes only the odd components (-90 degree, 90 degree) and the even components (+/-180 degree, 0 degree) from the 3D seismic data.

-----------------------------------------------------------------------------------------------------------------------------
**Scope of Program Field: Attribute Analysis - 3D Seismic Phase Decomposition**

_**Platform: MATLAB(R2020), Opendtect 6.4**_

Program file: "**Test_run_codes_pdoed3D_sample_3Dseismic_matfile.m**" provided script reads the sample seismic input data ("**sample_3Dseismic.mat**") and calls the phase decomposition function scripts except two files (**od_doprocess.m** & **od_getparameters.m**) and display a Inline results (./reference\images\Sample3Dseismic). User can change the input parameters values and display the result.

For Opendtect Volume Builder: Except "Test_run_codes_pdoed3D_sample_3Dseismic_matfile.m" files, all other MATLAB ".m" required to be complied via MATLAB C Compiler to generate the c shared dll file. How to setup the Environment, Opendtect-MATLAB link & to run the volume builder step from Opendtect, please refer "**./reference\Phase Decomposition Program Opendtect Link Setup.pdf**".

**Script does not support the batch processing in Opendtect and only valid for single 3D seismic data volume**

------------------------------------------------------------------------------------------
**References**:

* Case study: Phase-component amplitude variation with angle, GEOPHYSICS, VOL. 84, NO. 4 (JULY-AUGUST 2019); P. B285–B297, https://doi.org/10.1190/geo2018-0762.1

------------------------------------------------------------------------------------------
The usage of the program is limited to its scope and user has to ascertain the program output applicability to its scope of work, accuracy etc.

**Contact**: 
  
  Prashant Sinha [e-mail:sinha.pm@gmail.com]
