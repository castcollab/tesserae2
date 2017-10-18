import os
from setuptools import setup, find_packages
import json


def get_requirements_from_pipfile_lock(pipfile_lock=None):
    if pipfile_lock is None:
        pipfile_lock = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Pipfile.lock')
    lock_data = json.load(open(pipfile_lock))
    return [package_name for package_name in lock_data.get('default', {}).keys()]


packages = find_packages('.', exclude=['*.test', '*.test.*'])
pipfile_lock_requirements = get_requirements_from_pipfile_lock()

setup(
    name='pycortex',
    version='0.0.2',
    packages=packages,
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
    install_requires=pipfile_lock_requirements,
    python_requires=">=3.5",
)
