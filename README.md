# SAXS Z-score calculation


## System requirements:
1. Python 2.7
2. SciPy (tested with v1.0.0, any recent version should do)
3. NumPy (tested with v1.13.3, any recent version should do)
4. [IMP](https://integrativemodeling.org/) (optional, used to evaluate model-data fits, tested with v2.8.0)

## Installation guide
Place the files calc_z.py and curves.npz to a folder (optionally within your system PATH) and call as a python script. No specific configuration is required.

The online server which requires no installation is available [here](https://pharm.kuleuven.be/apps/biocryst/saxszs.php).

## Usage
python calc_z.py \<data\> [model]

Z0 score is calculated for an experimental scattering profile \<data\>.
If a PDB file [model] is provided, the script will additionally output chi2 score of the fit and the corresponding Z-score.


## Example

Example data can be found in the *test* folder. Run as "python calc_z.py original_0.dat original.pdb", it should take a couple of seconds on a modern PC. Expected outputs of the script:

#### Chi2 values

|SAXS profile   | original.pdb | perturbed.pdb |
|---------------|-------------:|--------------:|
|original_0.dat | 0.90         | 1.14          |
|original_1.dat | 0.89         | 1.93          |
|original_2.dat | 0.98         | 4.18          |
|original_3.dat | 1.04         | 10.06         |

#### Z-scores

|SAXS profile   | Z<sup>0</sup> | original.pdb | perturbed.pdb |
|---------------|--------------:|-------------:|--------------:|
|original_0.dat |   2.89        | 2.93         | 2.85          |
|original_1.dat |   3.34        | 3.37         | 3.14          |
|original_2.dat |   3.81        | 3.81         | 3.38          |
|original_3.dat |   4.15        | 4.14         | 3.47          |

The structural models and the SAXS profiles are based on the Tsirkone, V. G. et al.  J. Biol. Chem. 292, 9699–9710 (2017).
