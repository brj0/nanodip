# Verify running Python version (should be 3.7.5) and adjust jupyter notebook.
import IPython
import sys
from IPython.core.display import display, HTML
# set display witdth to 100%
display(HTML("<style>.container { width:100% !important; }</style>"))

sys.version

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
# end_internal_modules

"""
### ^^^ LIVE LOG ABOVE ^^^
All CherryPy access will be logged here, including live progress bars for
computationally intense analyses. Detailed access logging is turned off by
default (accessLogging is False), but can be turned on, e.g., for debugging,
in the configuration section at the beginning of this notebook. While it is not
required to have at look at these during normal operation, information
contained in the log may be helpful in troubleshooting. Line numbers in error
messages indicated here typically match those given in the respective Jupyter
Notebook cells.

To preserve these messages, halt the Python kernel, save and close the notebook
to send it for support. This makes sure that the code as well as the error
messages will be preserved.

To launch the user interface, wait until you see a pink log entry that the web
server has started, then navigate to http://localhost:8080.
"""
