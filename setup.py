from glob import glob
from io import open
from os.path import splitext, basename

from setuptools import setup, find_packages, Extension

# following src dir layout according to
# https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
version = "1.0.0"
setup(
    name="tesserae",
    version=version,
    description="Fast recombination-aware global and local alignment",
    author="Warren W. Kretzschmar, Kiran V Garimella",
    author_email="winni@warrenwk.com, kiran.garimella@gmail.com",
    license="Apache-2.0",
    long_description=open("README.rst").read(),
    install_requires="""
    numpy
    """.split(
        "\n"
    ),
    tests_require=["coverage", "pytest"],
    python_requires=">=3.6",
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    # entry_points={
    #     'console_scripts': ['cortexpy=cortexpy.__main__:main'],
    # },
    include_package_data=True,
)
