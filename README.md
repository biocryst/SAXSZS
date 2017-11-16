# SAXS Z-score calculation


## System requirements: 
1. Python 2.7
2. SciPy (tested with v1.0.0, any recent version should do)
3. NumPy (tested with v1.13.3, any recent version should do)
4. [IMP](https://integrativemodeling.org/) (optional, used to evaluate model-data fits, tested with v2.8.0)

## Installation guide:
Place the files calc_z.py and curves.npz to a folder (optionally within your system PATH) and call as a python script. No specific configuration is required.

The online server which requires no installation is available [here](https://pharm.kuleuven.be/apps/biocryst/saxszs.php).

## Usage:
python calc_z.py \<data\> [model]

Z0 score is calculated for an experimental scattering profile \<data\>.
If a PDB file [model] is provided, the script will additionally output chi2 score of the fit and the corresponding Z-score.


## Demo:

todo
