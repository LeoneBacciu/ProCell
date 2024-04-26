import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="procell",
    version="1.8.4",
    include_package_data=True,
    author="Marco S. Nobile",
    author_email="m.s.nobile@tue.nl",
    description="ProCell - cell proliferation framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aresio/ProCell",
    install_requires=['numpy', 'fst-pso', 'matplotlib', 'simpful', 'miniful', 'pyqt5', 'pandas'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
    ],
)