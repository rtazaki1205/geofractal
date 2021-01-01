# Introduction to geofractal

This package allows to compute average geometric cross sections of fractal dust aggregates 
by means of a statistical distribution model of monomers proposed in Tazaki (submitted to MNRAS).


# Terms of use

`geofractal` is distributed under the [MIE license](https://opensource.org/licenses/MIT) and can be used, changed
and redistributed freely. If you use this package to publish papers, please cite the following paper

> R. Tazaki
> *Analytical formulas of geometric cross sections of fractal aggregates*
> Submitted to MNRAS


# Examples 

## fortran

The input parameters can be set in `call.f90`.
The structure of fractal dust aggregates is specified by the fractal dimension `df`, and fractal prefactor `k0`, 
and number of monomers `N`. 
In addition to this, you also need to specify following three swithces 

- `iqapp` : Approximation for angular integration (=1 numerical solution, =2 approximate).
- `iqcon` : The numerical factor to connect small non-fractal cluster limit. 
- `iqcor` : The two-point correlation function of monomer distribution.

In default, `iqapp=2`,`iqcor=3`,`iqcon=2`.
To run the code, first, you make the codes
```
Make
```
This will create the executable file `results.x`, and then execute it by
```
./results.x
```
As a result, the output file `gratio.out` is created. 

## Python 

The python code assumes the recommended set of parameters: `iqapp=2`,`iqcor=3`,`iqcon=2`.
The input parameter can be set in `call.py`. To run the code, you need to perform

```
call call.py
```
As a result, the output file `gratio.out` is created. 



