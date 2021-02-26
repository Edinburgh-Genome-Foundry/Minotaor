from setuptools import setup, find_packages

version = {}
with open("minotaor/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="minotaor",
    version=version["__version__"],
    author="Peter Vegh",
    description="Amino acid sequence annotator",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/Minotaor",
    keywords="biology",
    packages=find_packages(exclude="docs"),
    install_requires=["pandas", "biopython"],
    include_package_data=True,
)
