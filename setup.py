import versioneer
from setuptools import setup


with open("README.md", "r") as src:
    LONG_DESCRIPTION = src.read()


setup(
    name="tsgxr",
    description="A tool for reading TSG spectra data into Xarray.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    url="https://github.com/CSIRO-GeoscienceAnalytics/tsg-xr",
    project_urls={
        "Code": "https://github.com/CSIRO-GeoscienceAnalytics/tsg-xr",
        "Issue tracker": "https://github.com/CSIRO-GeoscienceAnalytics/tsg-xr/issues",
    },
    author="Morgan Williams",
    author_email="morgan.williams@csiro.au",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    keywords=["TSG", "spectra"],
    packages=["tsgxr"],
    install_requires=["numpy", "xarray", "zarr", "pandas", "pytsg>=0.2.8", "typer"],
    extras_require={
        "plot": ["matplotlib"],
        "dev": ["pytest", "pytest-cov", "coverage", "versioneer", "black", "twine"],
    },
    tests_require=["pytest", "pytest-cov", "coverage"],
    test_suite="test",
    include_package_data=True,
    license="MIT",
    cmdclass=versioneer.get_cmdclass(),
    entry_points={
        "console_scripts": [
            "tsgxr = tsgxr.cli:app",
        ]
    },
)
