import setuptools

with open( "README.md" ) as fh:
    long_description = fh.read( )

setuptools.setup(
    name = "waafle",
    version = "0.1.0dev",
    author = "Eric Franzosa",
    author_email = "eric.franzosa@gmail.com",
    license = "MIT",
    description = "Identify LGT events from metagenomic contigs.",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://huttenhower.sph.harvard.edu/waafle",
    packages = setuptools.find_packages( ),
    classifiers = [
        "Programming Language :: Python",
    ],
    entry_points = {
        "console_scripts": [
            "waafle_genecaller = waafle.waafle_genecaller:main",
            "waafle_junctions = waafle.waafle_junctions:main",
            "waafle_orgscorer = waafle.waafle_orgscorer:main",
            "waafle_qc = waafle.waafle_qc:main",
            "waafle_search = waafle.waafle_search:main",
        ],
    },
    install_requires = [
        "numpy >= 1.13.0",
    ],
)
