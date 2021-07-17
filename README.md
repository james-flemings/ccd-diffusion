# Charged Coupled Device (CCD) Electron Diffusion

This repository contains a python script that runs a Monte Carlo Simulation of electrons diffusing in a CCD.

Two things you need to ensure so that the script works correctly

1. Install the neccessary packages. You'll need Numpy, Scipy, and Matplotlib. You can install them with the following command:

```bash
pip install -r requirements.txt
```

2. The mfp.txt is in the same directory as the mc.py

This Python script was written during my internship at Brookhaven National Laboratory in the Summer of 2020. I wrote an equivalent program in C++, but you'll need to install [ROOT](https://root.cern/install/) to run the program. By default, the Python script uses multiprocessing, so future work could be to make it optional (in case the user doesn't want to overwork his/her computer).
