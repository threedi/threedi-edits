from setuptools import setup

__version__ = "1.0"

long_description = "\n\n".join([open("README.rst").read(), open("CHANGES.rst").read()])

install_requires = ["scipy", "cached_property", "shapely=2.0a1", "threedigrid_builder"]

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
    description="Threedi Edits provides python tooling for threedi such as raster-conversion, alignment, fillers, checks and others",
    long_description=long_description,
    # Get strings from http://www.python.org/pypi?%3Aaction=list_classifiers
    classifiers=["Programming Language :: Python", "Framework :: Django"],
    keywords=[],
    author="Chris Kerklaan",
    author_email="chris.kerklaan@nelen-schuurmans.nl",
    url="https://github.com/threedi-edits",
    license="LGPL",
    packages=["threedi_edits"],
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requires,
    tests_require=tests_require,
    extras_require={"test": tests_require},
    entry_points={
        "console_scripts": ["run-threedi-edits = threedi_edits.scripts:main"]
    },
)
