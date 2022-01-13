# Installation

## Installing ShapePipe

The entire ShapePipe package, include dependencies, can be built as follows:

```bash
git clone https://drf-gitlab.cea.fr/cosmostat/ShapePipe
cd ShapePipe
./install_shapepipe
```

The `install_shapepipe` script will create the recommended Conda environment
along with all of the required core and module dependencies.

A list installation options can be seen using the `--help` option:

```bash
./install_shapepipe --help
```

## Installing the ShapePipe Library Only

The ShapePipe library, *i.e.* the core package not including module dependencies,
can be installed in the following ways.

### Using Git

The ShapePipe package can be downloaded from the repository
and built as follows:

```bash
git clone https://drf-gitlab.cea.fr/cosmostat/ShapePipe
cd ShapePipe
python setup.py install
```

This method is recommend for development.

### Using Pip

Alternatively, the ShapePipe library can be directly installed from the
repository as follows:

```bash
pip install git+ssh://git@drf-gitlab.cea.fr/cosmostat/ShapePipe.git
```

```{note}
Note, this method will not include any executable scripts or examples.
```

## Installing the Module Python Dependencies

Module Python dependencies can be installed in the following ways.

### Using Conda

The ShapePipe Conda environment can be built and activated by running:

```bash
conda env create -f environment.yml
source activate shapepipe
```

### Using Pip

Module Python dependencies can also be installed using `pip` as follows:

```bash
pip install -r requirements.txt
pip install -r requirements_git.txt
```
