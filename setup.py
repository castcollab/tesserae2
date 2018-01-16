import os
from setuptools import setup, find_packages
import json
from io import open

import cortexpy

version = cortexpy.__version__


def get_requirements_from_pipfile_lock(pipfile_lock=None):
    if pipfile_lock is None:
        pipfile_lock = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Pipfile.lock')
    lock_data = json.load(open(pipfile_lock))
    return [package_name for package_name in lock_data.get('default', {}).keys()]


packages = find_packages('.', exclude=['*.test', '*.test.*'])
pipfile_lock_requirements = get_requirements_from_pipfile_lock()

setup(
    name='cortexpy',
    version=version,
    description='The python sister project to CortexJDK',
    author='Warren W. Kretzschmar and Kiran Garimella',
    author_email='winni@warrenwk.com and kiran.garimella@gmail.com',
    maintainer='Warren W. Kretzschmar and Kiran Garimella',
    maintainer_email='winni@warrenwk.com and kiran.garimella@gmail.com',
    url='https://github.com/winni2k/cortexpy',
    download_url='https://github.com/winni2k/cortexpy/archive/{}.tar.gz'.format(version),
    license='Apache-2.0',
    long_description=open('README.rst').read(),
    install_requires=pipfile_lock_requirements,
    tests_require=['coverage', 'pytest'],
    python_requires=">=3.5",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    packages=packages,
    entry_points={
        'console_scripts': ['cortexpy=cortexpy.__main__:main_without_argv'],
    },
    include_package_data=True,
)
