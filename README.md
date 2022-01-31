# SSDFoam

## Compilation

```
./Allwmake
```

## Usage  

Include the compiled library to your ```controlDict```:

```
libs ( "libSSD.so" );
```

The scheme can be chosen in the ```fvSchemes``` dictionary as following:

```
div(phi,U)      GaussScaleSelective Gauss linear orthogonal 2.0 upwind blended 0.95
```

## References

<a id="1">[1]</a> 
Vuorinen, V., Larmi, M., Schlatter, P., Fuchs, L., & Boersma, B. J. (2012). A low-dissipative, scale-selective discretization scheme for the Navier–Stokes equations. Computers & Fluids, 70, 195-205, [10.1016/j.compfluid.2012.09.022](https://doi.org/10.1016/j.compfluid.2012.09.022).
<details>
<summary>BibTex</summary>
<p>
 
```
@article{vuorinen2012low,
  title={{A low-dissipative, scale-selective discretization scheme for the Navier--Stokes equations}},
  author={Vuorinen, Ville and Larmi, Martti and Schlatter, Philipp and Fuchs, Laszlo and Boersma, BJ},
  journal={Computers \& Fluids},
  volume={70},
  pages={195--205},
  year={2012},
  publisher={Elsevier}
}
```
 
</p>
</details>

<a id="2">[2]</a> 
Morev, I.: Implementation of scale-selective spatial discretization scheme in OpenFOAM. In Proceedings of CFD with OpenSource Software, 2021, Edited by Nilsson. H., [10.17196/OS_CFD#YEAR_2021](http://dx.doi.org/10.17196/OS_CFD#YEAR_2021).
<details>
<summary>BibTex</summary>
<p>
 
```
@inproceedings{morev2021implementation,
  author    = {Morev, Ilya},
  title     = {{Implementation of scale-selective spatial discretization scheme in OpenFOAM}},
  booktitle = {{Proceedings of CFD with OpenSource Software}},
  year      = {2021},
  publisher = {{Edited by Nilsson. H.}},
  doi={http://dx.doi.org/10.17196/OS_CFD#YEAR_2021}
}
```
 
</p>
</details>
