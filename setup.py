from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="plotly_phylotree",
    version="0.0.10",
    description="A package to create phylogenetic trees with Plotly",
    package_dir={"": "src"},
    packages=["phylotree"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/voelkerh/plotly.py.git",
    author="voelkerh",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=["plotly", "biopython", "numpy"],
    extras_require={
        "dev": ["pytest>=7.0", "twine>=4.0.2"],
    },
    python_requires=">=3.9",
)