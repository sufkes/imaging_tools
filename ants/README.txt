2020-03-06 - Steven Ufkes:
- This directory contains modified versions of ANTs scripts.

Changes:
antsRegistrationSyN-mutual_information.sh
* Use mutual information cost function instead of cross-correlation.

antsMultivariateTemplateConstruction2-modified_for_sickkids_hpc_pbs.sh
* Modified to work on the old PBS RIT-HPC cluster at SickKids.
* Added ANTSPATHALT which points to this directory.
* Use the alternate version of waitForPBSQJobs.pl (the one in this directory) everywhere.
* Add -l vmem=${MEMORY} to the qsub commands (only those called when the PBS qsub option is selected (option: -c 4)). This will hopefully prevent the vmem errors which have been occurring.
* Change args passed to waitForPBSQJobs.pl in registration step so that it polls qstat every 2 min instead of every 10 min (why not?)
* Probably does not work for multivariate template construction (only univariate).

waitForPBSQJobs-modified_for_sickkids_hpc_pbs.pl
* Modified to work on the old PBS RIT-HPC cluster at SickKids.
* Added logic so that jobs in state 'C' are considered complete.

antsMultivariateTemplateConstruction2-modified_for_ubc_sockeye.sh
* Various modifications to run on UBC ARC Sockeye cluster.
* Bug fixes to enable multivariate template construction.

waitForPBSQJobs-modified_for_ubc_sockeye.pl
* Various modifications to run on UBC ARC Sockeye cluster.
