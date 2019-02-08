from setuptools import setup, find_packages
from sphinx.setup_command import BuildDoc
cmdclass = {'build_sphinx': BuildDoc}

with open("README.md", "r") as handle:
    long_description = handle.read()

name = "genoml"
version = "0.0"
release = "0.0.1"

setup(
    name=name,
    version=release,
    author="Ross Altman",
    author_email="raltman@inari.com",
    description="Machine learning tools for genomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/knowbodynos/genoml",
    packages=find_packages(),
    cmdclass=cmdclass,
    # these are optional and override conf.py settings
    command_options={
        'build_sphinx': {
            'project': ('setup.py', name),
            'version': ('setup.py', version),
            'release': ('setup.py', release),
            'source_dir': ('setup.py', 'docs/source'),
            'build_dir': ('setup.py', 'docs/build')}},
    install_requires=[
        'six',
        'boto3',
        'botocore',
        'pandas',
        'pybedtools'
      ],
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
)