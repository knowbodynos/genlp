#!/bin/bash

# 1) Install sphinx (conda install -c anaconda sphinx, or pip install sphinx)
# 2) mkdir docs && cd docs
# 3) Run sphinx-quickstart, using
#    > Separate source and build directories (y/n) [n]: y
# 4) Add sys.path.insert(0, os.path.abspath('../../..')) to docs/source/conf.py
# 5) Run this script (bash ./build_docs.sh)

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

sphinx-apidoc -f -o ${DIR}/docs/source/ ${DIR}
python ${DIR}/setup.py build_sphinx