# Solvated OprD systems with different public small molecules

Note that these simulations are in explicit solvent with no membrane.

Files created by Camilo Velez-Vega, Entasis Therapeutics.

The pore axis is roughly aligned with the y-axis, with a corner of the water box at the origin.

Three residues in a plane to restrain: to prevent reorientation or measure height along pore:
* `GLY 342 CA`
* `ASP 97 CA`
* `SER 184 CA`

Center axis of pore (x,z) around `ALA 128 CA`.
Top of pore around (y-max) `SER 226 CA`.
Bottom of pore around (y-min) `VAL 88 CA`.

## Manifest

* `arg/` - arginine (compound 2 from [1])
* `glu/` - glutamate (not the same as compound 3 [glycine-glutamate] reported in [1])
* `his_posit` - histidine (positively charged, doubly protonated)
* `his_neut` - histidine (neutral, protonated at epsilon carbon)
* `imi/` - imipenem (compound 5 from [1])
* `mero/` - meropenem (compound 4 from [1])
* `comp7/` - compound 7 from [1]
* `comp8/` - compound 8 from [1]
* `comp7_nowat/` - compound 7 from [1] without solvent (for testing only)
* `comp8_nowat/` - compound 8 from [1] without solvent (for testing only)


## Running

Here's an example of how to run [`iapetus`](http://github.com/choderalab/iapetus) on a vacuum test system:
```bash
iapetus --gromacs comp7_nowat --ligseq 423 --output comp7_nowat.nc --n_steps_per_iteration 1250 --niterations 100 --verbose --testmode 
```
and on a real system:
```bash
iapetus --gromacs comp7 --ligseq 423 --output comp7.nc --n_steps_per_iteration 1250 --niterations 10000 --verbose
```


## FAQ

*Q:*  What do I do if I receive an error message like this?
```
Exception: CustomCentroidBondForce requires a device that supports 64 bit atomic operations
```
*A:* OpenMM is attempting to use a a device (such as an integrated laptop GPU) that doesn't support 64-bit atomic operations.
Either use a more powerful GPU, or if that is not possible, specify the CPU platform by adding the `--platform CPU` flag.

## References

[1] Isabella VM et al. Toward the Rational Design of Carbapenem Uptake in Pseudomonas aeruginosa. Chemistry & Biology 22:535, 2015. https://doi.org/10.1016/j.chembiol.2015.03.018

[2] Samanta S et al. Molecular basis of substrate translocation through the outer membrane channel OprD of Pseudomonas aeruginosa. PCCP 17:23867, 2015. http://doi.org/10.1039/C5CP02844B

[3] Iyer R et al. Whole-Cell-Based Assay To Evaluate Structure Permeation Relationships for Carbapenem Passage through the Pseudomonas aeruginosa Porin OprD. ACS Infect. Dis. 3:310, 2017. http://doi.org/10.1021/acsinfecdis.6b00197
