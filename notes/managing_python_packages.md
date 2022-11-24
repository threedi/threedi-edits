# Managing python packages on PyPi/test PyPi

## Building and uploading

1. Make sure your setup.py contains all the necessary information by running:

```python setup.py check```

If there is no output (or just ```running check```, all should be well.

2. Build a source distribution of the package:

From the root of the package (i.e. where setup.py is located), run:

```python setup.py sdist```

This creates a ```dist/``` directory (if it doesn't exist yet) and a compressed archive of the package.

3. Optionally you can create a wheel distribution (i.e. platform specific distribution, some parts are pre-built, requires `pip install wheel`)

```python setup.py bdist_wheel```

This creates a wheel in the ```dist/``` directory.

## Testing and publishing

### Uploading to testPyPi

```twine upload --skip-existing --repository-url https://test.pypi.org/legacy/ dist/*```

This uploads the newest build that isn't on testPyPi. You can then test the package installation by installing it
with the command listed for the project, something like:

```pip install -i https://test.pypi.org/simple/ [package_name==version]```
'''pip install --extra-index-url https://test.pypi.org/simple/ threedi_edits==0.13'''


This will not install dependencies specified in ```setup.py```, if you want this behaviour, add the following option to 
the command above:

```--extra-index-url https://pypi.org/simple``` before the package name and version.

### Uploading to PyPi

To upload a package to PyPi, simply use the following command:

```twine upload dist/[*]``` [*] -> choose the version and build you want to upload. e.g.: `twine upload dist/dist/hhnk_research_tools-0.2.0.tar.gz`.You can specify multiple files for 
the same upload


## Installing packages
`pip install hhnk-threedi-tests`

`pip install hhnk-research-tools`

```pip install --upgrade hhnk-threedi-tests --extra-index-url=https://test.pypi.org/simple/ ```



## Testing packges

python setup.py bdist_wheel

dist/hhnk_threedi_tests-0.1.9-py3-none-any.whl
pip install --upgrade dist/hhnk_threedi_tests-0.1.9-py3-none-any.whl
pip install --upgrade dist/hhnk_research_tools-0.2.1-py3-none-any.whl
restart qgis
