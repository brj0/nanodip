#!/usr/bin/env python

"""
## NanoDiP all-in-one Jupyter Notebook
*J. Hench, S. Frank, C. Hultschig and J. Brugger, Neuropathology, IfP Basel,
2021*

This software is provided free of charge and warranty; by using it you agree
to do this on your own risk. The authors shall not be held liable for any
damage caused by this software. We have assembled this and tested it to the
best of our knowledge.

The purpose of NanoDiP (Nanopore Digital Pathology) is to compare low-coverage
Nanopore sequencing data from natively extracted DNA sequencing runs against
a flexibly adaptable collection of 450K/850K Illumina Infinium Methylation
array data. These data have to be preprocessed into binary beta value files;
this operation is performed in R (uses minfi to read raw array data) and
outputs bindary float files (one per dataset). These beta values files (e.g.,
204949770141_R03C01_betas_filtered.bin) are named according to the array ID
(Sentrix ID) followed by the suffix. A collection of betas_filtered.bin files
can be provided in a static manner and XLSX (Microsoft Excel) tables can be
used to select a subset thereof alongside a user-defined annotation. The
corresponding datasets will be loaded into memory and then serve as the
reference cohort to which the Nanopore data are compared by dimension reduction
(UMAP). This comparison is optimized for speed and low resource consumption so
that it can run on the computer that operates the sequencer. The sequencing run
is initiated through the MinKNOW API by this application. Basecalling and
methylation calling occur as background tasks outside this Jupyter Notebook.
User interaction occurs through a web interface based on CherryPy which has
been tested on Chromium web browser. It is advisable to run it locally, there
are no measures to secure the generated website.

In order to use this application properly please make sure to be somewhat
familiar with Jupyter Notebook. To run the software, press the button called
*restart the kernel, re-run the whole notebook (with dialog)* and confirm
execution. Then, in Chromium Browser, navigate to http://localhost:8080/ and
preferably bookmark this location for convenience. In case of errors, you may
just again click the same button *restart the kernel, re-run the whole notebook
(with dialog)*.
___

### Technical Details

* Tested with Python 3.7.5; 3.8.8 fails to load minknow_api in jupyter
  notebook.
* Verified to run on Ubuntu 18.04/Jetpack on ARMv8 and x86_64 CPUs; not
  tested on Windows and Mac OS. The latter two platforms are unsupported, we
  do not intend to support them.
* **CAUTION**: Requires a *patched* version of minknow api, file
  `[VENV]/lib/python3.7/site-packages/minknow_api/tools/protocols.py`.
  Without the patch, the generated fast5 sequencing data will be unreadable
  with f5c or nanopolish (wrong compression algorithm, which is the default in
  the MinKNOW backend).
"""

# Verify running Python version (should be 3.7.5) and adjust jupyter notebook.
import IPython
import os
from IPython.core.display import display, HTML
# set display witdth to 100%
display(HTML("<style>.container { width:100% !important; }</style>"))
os.system('python --version')

"""
## Multithreading Options
Depending on the number of parallel threads/cores of the underlying hardware,
threading options for multithreaded modules need to be set as
environment-specific parameters. One way to do so is through the *os* module.
"""

# execution-wide multithreading options, set according to your hardware. Jetson
# AGX: suggest "2" needs to be set before importing other modules that query
# these parameters
import os
os.environ["NUMBA_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"

"""
## Modules
This section imports the required modules that should have been installed via
pip. Other package managers have not been tested. To install packages, use the
setup script provided with this software or, alternatively, install them one
by one, ideally in a virtual python environment. Note that the MinKNOW API
requires manual patching after installation with pip.
"""

# python_modules_to_import
# start_external_modules
# end_external_modules

# start_internal_modules
from webui import (
    start_webserver,
)
# end_internal_modules

if __name__ == "__main__":
    start_webserver()

"""
### ^^^ LIVE LOG ABOVE ^^^
All CherryPy access will be logged here, including live progress bars for
computationally intense analyses. Detailed access logging is turned off by
default (accessLogging is False), but can be turned on,e.g., for debugging,
in the configuration section at the beginning of this notebook. While it is not
required to have at look at these during normal operation, information
contained in the log may be helpful in troubleshooting. Line numbers in error
messages indicated here typically match those given in the respective Jupyter
Notebook cells.

To preseve these messages, halt the Python kernel, save and close the notebook
to send it for support. This makes sure that the code as well as the error
messages will be preserved.

To launch the user interface, wait until you see a pink log entry that the web
server has started, then navigate to http://localhost:8080.
"""
