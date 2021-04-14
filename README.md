# Introduction to geofractal

This package allows to compute average geometric cross sections of fractal dust aggregates 
by means of a statistical distribution model of monomers proposed in Tazaki (submitted to MNRAS).  

# Terms of use

`geofractal` is distributed under the [MITlicense](https://opensource.org/licenses/MIT) and can be used, changed
and redistributed freely. If you use this package to publish papers, please cite the following paper

> R. Tazaki  
> * Analytic expressions for geometric cross-sections of fractal dust aggregates*  
> Accepted for publication in MNRAS  


# Examples 

## fortran 

The input parameters can be set in `call.f90`.  
The user can specify following three input parameters in `call.f90`
- `df` : Fractal dimension (1 ≦ df ≦ 3)
- `k0` : Fractal prefactor
- `PN` : Number of monomers (1 ≦ PN)

In addition, the user also needs to specify following three options

- Angular integration in the calculation of the mean overlapping efficiency  
	`iqapp=1` : use numerical integration  
	`iqapp=3` : use approximate analytical solution   
- Small non-fractal cluster limit  
	`iqcon=1` : without small cluster limit  
	`iqcon=2` : with small cluster limit  
- The two-point correlation function of fractal aggregates  
 	`iqcor=1` : The Gaussian cut-off model  
	`iqcor=2` : The exponential cut-off model  
	`iqcor=3` : The fractal dimension cut-off model  
	
I recommend following set of options: `iqapp=3`,`iqcor=3`,`iqcon=2` (default).  

To run the code, first of all, compile the codes via
```
make
```
This will create an executable file `results.x`. Then, perform
```
./results.x
```
As a result, the output file `gratio.out` is created. 

## python 

The python package of `geofractal` is a shortened version of the fortran package as it runs only with the recommended options: `iqapp=3`,`iqcon=2`.

Similar to the fortran package, the user can specify the input parameters in `call.py`
- `df` : Fractal dimension (1 ≦ df ≦ 3)
- `k0` : Fractal prefactor
- `PN` : Number of monomers (1 ≦ PN)

In addition, the use also needs to opt a model of two-point correlation function:  
- `cormodel='EXPNL'` :  The exponential cut-off model  
- `cormodel='GAUSS'` :  The Gaussian cut-off model  
- `cormodel='FLDIM'` :  The fractal dimension cut-off model (default)

To run the code, simply perform
```
python call.py
```
As a result, the output file `gratio.out` is created. 


# History

*Version 0.1*
- Pre-release
