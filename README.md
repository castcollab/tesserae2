[![Build Status](https://travis-ci.com/winni2k/pycortex.svg?token=K7dhHdBzXsubBntxA949&branch=master)](https://travis-ci.com/winni2k/pycortex)

# pycortex
Python implementation of [CortexJDK](https://github.com/mcveanlab/CortexJDK)

## Development

Install pipenv:
```bash
pip install pipenv
```

Install dependencies
```bash
pipenv install
```

Activate environment
```bash
pipenv shell
```

## Todo

#### Finish Cortex Graph
- [X] Implement straight up iteration over .ctx file
- [X] Implement log(n) random access to sorted .ctx file
- [X] Implement CortexJDK Print method

#### Provide Cortex Graph adapter for use with networkx
