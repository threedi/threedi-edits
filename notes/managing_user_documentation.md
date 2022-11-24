# Managing user documentation

To manage the user documentation, activate the corresponding .yml file. If you get permission errors, you can 
manually install the dependencies with the ```--user``` option.

Then navigate into de 
```docs``` directory where the user manual is stored.

To rebuild the docs, run:

```make clean``` then ```make html``` (other builds are possible, however, there are hardcoded html 
links in the source files which won't translate (sorry! :|))

If all is well, ```readthedocs``` automatically tracks the git repository where the files are stored, 
so once changes are pushed they should be reflected. It is probably wise to check when the last 
build was to make sure.

