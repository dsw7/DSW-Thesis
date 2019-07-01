## Met-aromatic low level routine
---  
A low level object oriented package for running the Met-aromatic algorithm. See DSW thesis for a theoretical description. Met-aromatic data is first collected using ```runner.py``` which loads data into a MongoDB database. The script ```analyze.py``` runs a set of queries on the MongoDB database and returns meaningful data which can be redirected for storage. The script ```heatmap.py``` creates a heatmap representation of bridging interactions broken down by EC classifiers. This was a reviewer request.

### Usage: runner.py
---
To run the program:
```
$ python runner.py <args>
```
The program requires, at bare minimum, either a valid PDB code or a path to a text file containing a list of PDB codes to analyze. Here the arbitrary PDB code 1rcy is analyzed:
```
$ python runner.py --code 1rcy
```
A batch job can be performed as follows:
```
$ python runner.py --batch /path/to/foo.txt
```
The PDB codes in the batch job text file should be separated by newline characters. *NOTE:* Both code and batch parameters cannot be passed simultaneously. Next come the Met-aromatic algorithm constraints:
```
$ python runner.py --code 1rcy --cutoff 4.9 --angle 90.0 --model cp
```
Here the cutoff has been set to 4.9 Angstroms (the max norm of vector *v*) and the maximum angle of either Met-theta or Met-phi cannot exceed 90.0 degrees. The model used to interpolate lone pair positions is cp or Cross Product. These parameters do not have to be passed. Default values are used if these values are not specified. Defaults can be obtained by reading:
```
$ python runner.py --help
```
Console output is normally suppressed. Suppression can be lifted by passing the verbose parameter:
```
$ python runner.py --code 1rcy --verbose
```
There are several options available for working with output data. Data can either be exported to a .csv file or a MongoDB database. Export to a MongoDB database is recommended. Data cannot be exported to both a .csv file and a MongoDB database simultaneously. To save to a .csv file:
```
$ python runner.py --code 1rcy --export-csv /path/to/output.csv
```
Or a MongoDB database:
```
$ python runner.py --code 1rcy --export-mongo
```
MongoDB export can be modified as follows:
```
$ python runner.py --code 1rcy --export-mongo --mongoport 27017 --mongohost localhost --database my_database --collection my_collection
```
Default MongoDB parameters are passed if no export parameters are specified. No data is saved if no export parameter is passed. As always, defaults can be obtained using:
```
$ python runner.py --help
```

### Usage: analyze.py
---
To run the script:
```
$ python analyze.py
```
Script interpretation will be terminated if data was not previously collected using ```runner.py```. Otherwise, data will be printed to the console. Output can be redirected to a file for storage:
```
$ python analyze.py > /path/to/results/results.txt
```

### Usage: heatmap.py
---
There's really little to this script. I literally copy pasted data from ```analyze.py``` into this script.
```
$ python heatmap.py
```
Which yields:
<img src="https://github.com/dsw7/DSW-Thesis/blob/master/thesis_Met-aromatic_workdir_LOWLEVEL/2_bridges.png", align="center">
