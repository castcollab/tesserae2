import os
import sys
import json
from io import open
from setuptools import setup, find_packages
from Cython.Build import cythonize

sys.path.append(os.path.join(__file__, 'cortexpy'))
from cortexpy import __version__

version = __version__


def get_requirements_from_pipfile_lock(pipfile_lock=None):
    if pipfile_lock is None:
        pipfile_lock = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Pipfile.lock')
    lock_data = json.load(open(pipfile_lock))
    return [package_name for package_name in lock_data.get('default', {}).keys()]


packages = find_packages('.')
pipfile_lock_requirements = get_requirements_from_pipfile_lock()

setup(
    name='cortexpy',
    version=version,
    description='The python sister project to CortexJDK',
    author='Warren W. Kretzschmar, Kiran V Garimella',
    author_email='winni@warrenwk.com, kiran.garimella@gmail.com',
    url='https://github.com/winni2k/cortexpy',
    download_url='https://github.com/winni2k/cortexpy/archive/{}.tar.gz'.format(version),
    license='Apache-2.0',
    long_description=open('README.rst').read(),
    install_requires=pipfile_lock_requirements,
    tests_require=['coverage', 'pytest'],
    python_requires=">=3.6",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    packages=packages,
    entry_points={
        'console_scripts': ['cortexpy=cortexpy.__main__:main_without_argv'],
    },
    include_package_data=True,
    ext_modules=cythonize("cortexpy/graph/parser/kmer_ext.pyx")#, annotate=True)
)
