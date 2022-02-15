# Information regarding the rhoLENT data

* The parameter variation is prepared with PyFoam:  

    * .parameter files define the Cartesian products of parameters.
    * .variation files connect the simulation-ID in the format "%05d" with a parameter vector.
    * .json files agglomerate parameter variation data, they are read with `pandas.read_json(filename, orient='table)` 

* The relevant parameter variation files are the *.html exports of Jupyter Notebooks and the *.ipynb notebook files themselves. 
    * The python modules required for visualization are available in the rhoLENT-data folder so the Jupyter notebooks are self-sustained.
