# Constraining national NOx emission totals and their trends using satellite data

Measuring NOx concentration from space is relatively easy. The pollutant's lifetime in the atmosphere is short, non-anthropogenic sources are few and absorption patterns are distinct. Thus, NOx is the perfect candidate for a top-down emission inventory approach. Moreover, satellite data on atmospheric NOx concentration is available from about 1995 onwards, also allowing for trend analysis and comparison to the official bottom-up emission inventory figures. For this project, past and current satellite data will be examined and used to constrain national NOx emission estimates both regarding their total values and their trends. The methodology used should be applicable for other pollutants, such as NH3 and SO2. To this end, particular focus should be given to emission estimates for area sources as they are commonly the sources associated with the highest uncertainties and statistical difficulties.

# Tool

We aim at creating and then freely offering a software tool here. The program allows for emissions inventory estimates to be compared with independently derived calculations based on satellite data. Please refer to the repository's [wiki](https://github.com/UBA-DE-Emissionsituation/space-emissions/wiki) for all the details.

# Setting up a working python environment

We use specific python packages within these tools. If the tool does not run off the shelf and errors suggest this has to do with missing/failing packages, making a dedicated conda environment can be the solution. To enable usage of the same environment in which the tools have been tested a YAML file (environment.yml) is inculded in this project. Generating a conda environment, called geo_env, with this file is done using "conda env create --file environment.yml". Don't forget activating the enviroment before using the tool through "conda activate geo_env".
