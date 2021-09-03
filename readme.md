# Receptor Builder

### Env
```shell
conda create -n receptorbuilder python=3.7 -y 
conda activate receptorbuilder
conda install -c openeye openeye-toolkits=2020 -y 
conda install -c conda-forge fpocket -y 
conda install numpy 
```

### Usage

`python main.py mypdbfile.pdb`

This will create a folder called `mypdbfile_out`, and inside it you
will find the results of FPocket along with four .oeb files. Only the four
biggest and high scoring pockets are built into oeb files since most of the others are too small.

If you want to build all the pockets, remove lines 133 and 134 in main.