NanoDiP
===============================================================================


<h1 align="center">
<img src="/nanodip/static/img/logo.svg" width="300">
</h1><br>

This project is a fork of https://github.com/neuropathbasel/nanodip and still **buggy and under construction**.


Disclaimer
-------------------------------------------------------------------------------
Implementation of this software in a diagnostic settings occurs in the sole responsibility of the treating physician. Usage of this software occurs at the risk of the user. The authors may not be held liable for any damage (including hardware) this software might cause.


Installation
-------------------------------------------------------------------------------
The intended location for the installation is `/applications/nanodip`. For a installation within a virtual environment you can use:

```sh
mkdir -p /applications/nanodip_env
virtualenv --python 3.7 /applications/nanodip_env
git clone https://github.com/brj0/nanodip /applications/nanodip
cd /applications/nanodip
source /applications/nanodip_env/bin/activate
pip install -U setuptools
pip install -e .
```


NanoDiP can be started with `python /applications/nanodip/nanodip/main.py` or simply with the console command `nanodip`.
