[![Build Status](https://travis-ci.org/winni2k/cortexpy.svg?branch=master)](https://travis-ci.org/winni2k/cortexpy)
[![Coverage Status](https://coveralls.io/repos/github/winni2k/cortexpy/badge.svg?branch=admin_category)](https://coveralls.io/github/winni2k/cortexpy?branch=admin_category)
[![Code Health](https://landscape.io/github/winni2k/cortexpy/master/landscape.svg?style=flat)](https://landscape.io/github/winni2k/cortexpy/master)

# Cortexpy

Python implementation of [CortexJDK](https://github.com/mcveanlab/CortexJDK)

Cortexpy supports python versions 3.5 and 3.6.

## Development

Set up the project

```bash
make init
```

See `cortexpy/__main__.py` for example uses of cortexpy.

### Publishing to Pypi

```bash
# assuming version 0.0.1

# Step 1: Set __version__ in cortexpy/__init__.py

# Step 2: Make version commit
git commit -am'Bump to version 0.0.1'

# Step 3: tag commit
git tag 0.0.1

# Step 4: publish to pypi
make deploy
```
