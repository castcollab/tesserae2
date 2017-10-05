from setuptools import setup, find_packages

setup(
    name='pycortex',
    version='0.0.1',
    packages=find_packages('pycortex'),
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
    python_requires=">=3.5",
)
