# greta
GRoupEd sTar formAtion: converting sink particles to stars that are sampled from observational initial mass function (IMF).

Paper: Liow K. Y., Rieder S., Dobbs C. L., Jaffa S. E., 2022, MNRAS, 510, 2657 (ADS link of the paper [here](https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.2657L/abstract)).

Grouped star formation is developed to form stars dynamically in parsec-scale simulations, such that that the sink particles are first grouped together, then the total group masses are used to sample the IMF and form stars. The current version is not designed to be used to generate initial conditions.

`greta` is built upon the [AMUSE framework](https://github.com/amusecode/amuse).

# Installation

First, install `AMUSE`:

```bash
pip install amuse-framework
pip install amuse
```

Then, install my version of `AMUSE MASC` (c.f. Steven Rieder's [MASC](https://github.com/rieder/MASC)) and create a symbolic link:

```bash
git clone https://github.com/kyliow/MASC.git masc
cd {amuse directory}/src/amuse/ext
ln -s ../../../../masc/src/amuse/ext/masc
cd ../../../../
```  

Lastly, install `greta`:

```bash
git clone https://github.com/kyliow/greta.git greta
```

# How to use

The command line arguments that can be parsed to `greta.py` are:
- `-i`: Input AMUSE sink particle set file name.
- `-o`: Output AMUSE star particle set file name. Default is `stars.amuse`.
- `-d`: Grouping distance parameter in pc. Default is 0.
- `-v`: Grouping speed parameter in km/s. Default is 0.
- `-t`: Grouping age parameter in Myr. Default is 0.
- `--lower-limit`: Lower limit of stellar mass in MSun. Default is 0.5.
- `--upper-limit`: Upper limit of stellar mass in MSun. Default is 100.
- `-r`: Random seed. Default is `None`.

To 'turn off' any of the three grouping parameters, simply set them to unrealistically high values. See Liow et al. (2022) for more information.

An example command to run `greta.py` with grouping parameters 1 pc, 3 km/s and 1 Myr:

```bash
python greta.py -i sinkfilename.amuse -o output.amuse -d 1 -v 3 -t 1
```
