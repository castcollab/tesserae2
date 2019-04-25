from glob import glob
from io import open
from os.path import splitext, basename

from setuptools import setup, find_packages, Extension

try:
    from Cython.Build import cythonize
except ImportError:
    USE_CYTHON = False
else:
    USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [Extension(
    "cortexpy.graph.parser.kmer_ext",
    ["src/cortexpy/graph/parser/kmer_ext" + ext],
    extra_compile_args=["-std=c++11"],
    extra_link_args=["-std=c++11"]
)]

if USE_CYTHON:
    extensions = cythonize(extensions)  # , annotate=True)

version = '0.46.4'
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
    install_requires="""
    attrs
    biopython
    numpy
    networkx==2.1
    schema
    delegation
    msgpack
    pyyaml>=4.2b1
    """.split('\n'),
    tests_require=['coverage', 'pytest'],
    python_requires=">=3.6",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    ext_modules=extensions,
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
    entry_points={
        'console_scripts': ['cortexpy=cortexpy.__main__:main'],
    },
    include_package_data=True
)
