### Met-aromatic LOW LEVEL routine
A low level object oriented package for running the Met-aromatic algorithm. See DSW thesis for a theoretical description.  

To run the program:
```
$ python runner.py <args>
```
The program requires, at bare minimum, either a valid PDB code or a path to a text file containing a list of PDB codes to analyze:
```
$ python runner.py --code 1rcy
```
Or:
```
$ python runner.py --batch /path/to/foo.txt
```
The PDB codes in the text file should be separated by newline characters. *NOTE:* Both parameters cannot be passed simultaneously. Next come the Met-aromatic algorithm constraints:
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





