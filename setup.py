from setuptools import setup, find_packages


__version__ = "1.2.2"

long_description = "\n\n".join([open("README.rst").read(), open("CHANGES.rst").read()])

install_requires = [
    "scipy",
    "cached_property",
    "shapely==2.0a1",
    "threedigrid_builder",
    "threedi_modelchecker",
]

tests_require = [
    "pytest",
    "mock",
    "pytest-cov",
    "pytest-flakes",
    "pytest-black",
    "pytest-benchmark",
]

setup(
    name="threedi-edits",
    version=__version__,
    description="An experimental pythonic 3Di schematisation api. Using this api, we can access, alter en write a 3Di database within python. Within the package gis tools are provided as well. ",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    # Get strings from http://www.python.org/pypi?%3Aaction=list_classifiers
    classifiers=["Programming Language :: Python", "Framework :: Django"],
    keywords=["3Di", "GDAL", "api"],
    author="Chris Kerklaan",
    author_email="chris-kerklaan@live.nl",
    url="https://github.com/threedi-edits",
    license="LGPL",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requires,
    tests_require=tests_require,
    python_requires="<3.11",
    extras_require={"test": tests_require},
    entry_points={
        "console_scripts": ["run-threedi-edits = threedi_edits.scripts:main"]
    },
)
