import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="genoml",
    version="0.0.1",
    author="Ross Altman",
    author_email="raltman@inari.com",
    description="Machine learning tools for genomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/knowbodynos/genoml",
    packages=setuptools.find_packages(),
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
)