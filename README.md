# imaging_tools
Scripts for various image processing tasks

## How to run scripts (for Miller lab members)
Members of the Miller lab with an RIT-HPC account should be able to run these scripts on HPC in interactive sessions or batch scripts, after running the following command:

```source /hpf/largeprojects/smiller/tools/python_environments/activate_py3env.sh```

This will load the required HPC modules, and will activate a Python environment with the required Python packages installed.

### `dicom/` - Tools for processing DICOM files
### `dti/` - Tools for processing and analyzing diffusion tensor imaging data
### `fmri/` - Tools for processing and analyzing functional MRI data
### `minc/` - Tools for processing MINC files
### `misc/` - Tools for miscellaneous image processing tasks 

## Requirements
The following software is required. All packages are already installed on HPC and can be activated as described above.
* Python 3.9
* Python 2.7 (for some old scripts).
* Python packages:
  * numpy
  * scipy
  * matplotlib
  * pandas
  * pydicom
  * pyminc
  * nibabel
  * natsort
  * bctpy
  * sklearn
  * nilearn
* DCMTK (https://dicom.offis.de/index.php.en)
* FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
* MINC Toolkit (https://bic-mni.github.io/)