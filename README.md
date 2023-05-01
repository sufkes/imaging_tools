# imaging_tools
Scripts for various image processing tasks

## How to run scripts (for Miller lab members)
Members of the Miller lab with an RIT-HPC account should be able to run these scripts on HPC in interactive sessions or batch scripts, after running the following command:

```source /hpf/largeprojects/smiller/tools/python_environments/activate_py3env.sh```

This will load the required HPC modules, activate a Python environment with the required packages installed, and add the script directories to your `$PATH` variable so that they can be run by name.

## Available tools
Scripts that do not have an entry in this list may not work yet.
### `dicom/` - Tools for processing DICOM files
* `processDicom.py` - Sort DICOM files into directories.
  * Run `processDicom.py -h` for usage.
* `anonDicom.py` - Deidentify DICOM files.
  * Run `anonDicom.py -h` for usage.
### `dti/` - Tools for processing and analyzing diffusion tensor imaging data
### `fmri/` - Tools for processing and analyzing functional MRI data
### `misc/` - Tools for miscellaneous image processing tasks
* `makeThumbs.py` - Generate thumbnails of 3D images with overlays for rapid viewing.
  * Run `makeThumbs.py -h` for usage.
* `padNifti.py` - Add zero padding slices/volumes to NIFTI image. Padding can be added to specific sides of specific dimensions.
  * Run `padNifti.py -h` for usage.
* `pc2da.py` - Convert Canadian postal codes to dissemination area unique identifiers (DAuid).
  * Run `pc2da.py -h` for usage.
  * Input a postal code conversion file (https://mdl.library.utoronto.ca/collections/numeric-data/census-canada/postal-code-conversion-file) and return a spreadsheet containing all postal codes and corresponding DAuids.
  * Optionally, input an additional spreadsheet with single column of postal codes, and convert only those postal codes listed.
* `roiStats.py` - Report statistics on size and extent of region of interest clusters in 3D images.
  * Run `roiStats.py -h` for usage.
  * For examples, see `/hpf/largeprojects/smiller/examples/lesion_extent/README.txt` on HPC.
### `old_scripts/` - Old scripts that are kept for reference.

## Requirements
The following software is required. All packages are already installed on HPC and can be activated as described above.
* Python 3.9 or newer.
* Python 2.7 (for some of the older scripts).
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
  * scikit-image
* DCMTK (https://dicom.offis.de/index.php.en)
* FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
* MINC Toolkit (https://bic-mni.github.io)
* ANTs (https://github.com/ANTsX/ANTs)

### Example creation of Python 2.7 and 3.11 environments with all module requirements using Conda for environment creation and Pip for package installation. 
```
# Create Python 3.11 environment. 
conda create --prefix <path to conda environment directory>/py3 python=3.11
source activate <path to conda environment directory>/py3
pip install numpy scipy matplotlib pandas nibabel pyminc Pillow opencv-python pydicom bctpy scikit-learn scikit-image torch torchvision torchaudio PyCap

# Create Python 2.7 environment
conda create --prefix <path to conda environment directory>/py2 python=2.7
source activate <path to conda environment directory>/py2
pip install numpy scipy matplotlib pandas nibabel pyminc Pillow pydicom pyyaml pycap==1.0.2
```