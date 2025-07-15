# Exoclimes VII 3D clouds practical
 
Python based practical for Exoclimes VII cloud modelling review talk.
This repository acts as a supplemental to the presentation (included here too).

Each example code has extended notes and line-by-line comments that detail each part of each method. 
Suggestions, ideas and improvements that help the learning/user experience are most welcome!

Each example uses a 1D approach to demonstrate the salient features of each approach.
It is relatively simple to design 3D versions and couple to hydrodynamic models (e.g. GCMs).

* Example 1: Ackerman & Marley (2001) model
* Example 2: Tracer saturation adustment model
* Example 3: Two-moment monodisperse mass moment mixed material model

NOTE: The best way to generate a non-irradiated profile (e.g. brown dwarf) is set the zenith angle to a small value (e.g. mu_z = 1e-6) and Tirr to a low value (e.g. Tirr = 1), rather than set Tirr = 0. 
Then set the Tint to the required Teff. 
This is just due to the (probably bad) way the semi-grey and picket-fence T-p profile calculation was coded.

### Running the code

To run the code for each example, for example, example 1.py, enter the following in the terminal:

> python example_1.py

You will need to install the following packages for the examples to work, (e.g. using conda or pip):

+ pyyaml
+ numpy
+ scipy
+ matplotlib
+ seaborn
+ numba (for @jit constructs)

This should generate matplotlib plot popups, as well as png files in the directory with the results.

Parameters for each example can be modified in the parameters.yaml file - feel free to modify this to try out different scenarios.
Various condensation species can be set, see the `vapour_pressure' function in the atm_module.py file for a complete list.
Then set the appropriate paramaters in the parameter.yaml file for that species.

WARNING: examples 2 and 3 will take longer to calculate than 1, in the small amount of integration time used as the default here, you're unlikely to achieve a converged solution.


## Key papers for different methods (non-exhaustive)

### Tsuji

* [Tsuji (2002)](https://ui.adsabs.harvard.edu/abs/2002ApJ...575..264T/abstract)
* [Tsuji (2005)](https://ui.adsabs.harvard.edu/abs/2005ApJ...621.1033T/abstract)

### Rossow + Allard (BT-Settl)

* [Rossow (1978)](https://ui.adsabs.harvard.edu/abs/1978Icar...36....1R/abstract)
* [Allard et al. (2001)](https://ui.adsabs.harvard.edu/abs/2003IAUS..211..325A/abstract)
* [Allard et al. (2007)](https://ui.adsabs.harvard.edu/abs/2007A%26A...474L..21A/abstract)

### Ackerman & Marley (EddySed & VIRGA)

* [Ackerman & Marley (2001)](https://ui.adsabs.harvard.edu/abs/2001ApJ...556..872A/abstract)
* [Rooney et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...33R/abstract)

### Ormel & Min (ARCiS)

* [Ormel & Min (2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...622A.121O/abstract)
* [Huang et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024A%26A...691A.291H/abstract)

### Helling & Woitke (DRIFT & DIFFU-DRIFT)

* [Woitke & Helling (2003)](https://ui.adsabs.harvard.edu/abs/2003A%26A...399..297W/abstract)
* [Woitke & Helling (2004)](https://ui.adsabs.harvard.edu/abs/2004A%26A...414..335W/abstract)
* [Helling & Woitke (2006)](https://ui.adsabs.harvard.edu/abs/2006A%26A...455..325H/abstract)
* [Helling et al. (2008)](https://ui.adsabs.harvard.edu/abs/2008A%26A...485..547H/abstract)
* [Woitke et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020A%26A...634A..23W/abstract)

### Tracer saturation adjustment (Used in various GCMs)

* [Tan & Showman (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...874..111T/abstract)
* [Tan & Showman (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.2198T/abstract)
* [Komacek et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...934...79K/abstract)
* [Lee et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529.2686L/abstract)

### Mass moment method (mini-cloud)

* [Ohno & Okuzumi (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...835..261O/abstract)
* [Lee & Ohno (2025)](https://ui.adsabs.harvard.edu/abs/2025A%26A...695A.111L/abstract)
* [Lee (2025)](https://ui.adsabs.harvard.edu/abs/2025arXiv250310309L/abstract)
* [Lee & Ohno (A&A submitted)]()

### Bin methods (CARMA)

* [Powell et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...860...18P/abstract)
* [Gao et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...855...86G/abstract)
* [Gao et al. (2018b)](https://ui.adsabs.harvard.edu/abs/2018ApJ...863..165G/abstract)
* [Powell et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...887..170P/abstract)
* [Powell et al. (2022)](https://ui.adsabs.harvard.edu/abs/2024ApJ...969....5P/abstract)


## DISCLAIMER:

These are codes for learning basics and toy modelling only, and are NOT complete with regards to physics/chemistry in-depth details (which matter a lot), may contain bugs, and are not as tested as the original implementations. They are also not numerically optimised (on purpose for readability).

So do NOT use these for end product science, consult an expert in the field for example (non-exhaustive list) A&M -> VIRGA [Natasha Batalha, Caroline Morley, Mark Marley], tracer sat adj -> [Xianyu Tan, Tad Komacek, Elspeth Lee], two mass moment microphysics -> [Kazumasa Ohno, Elspeth Lee] if you want to go further with any of the methods for `real' science.

For bin models (not used here) consult a CARMA expert e.g. Diana Powell, Peter Gao. \
For interest in the DRIFT moment methodology consult Christiane Helling or Peter Woitke. \
For interest in advanced size distribution dependent mass moment methods consult Elspeth Lee or Kazumasa Ohno.
