"""Setup configuration for uvlf-hod package."""

from setuptools import setup, find_packages
import os

# Read version
version = {}
with open("uvlf_hod/__version__.py") as f:
    exec(f.read(), version)

# Read README for long description
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="uvlf-hod",
    version=version["__version__"],
    author="Marko Shuntov",
    author_email="marko.shuntov@nbi.ku.dk",
    description="UV Luminosity Function and Halo Occupation Distribution modeling for high-redshift galaxies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mshuntov/uvlf-hod",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=3.0",
            "black>=22.0",
            "flake8>=4.0",
            "sphinx>=4.0",
        ],
    },
    keywords="astronomy cosmology galaxies luminosity-function halo-model",
)
