#!/bin/bash

# 1) Install sphinx (conda install -c anaconda sphinx, or pip install sphinx)
# 2) mkdir docs && cd docs
# 3) Run sphinx-quickstart
# 4) Run this script (bash ./build_docs.sh)

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

sphinx-apidoc -f -o ${DIR}/docs/source/ ${DIR}
python ${DIR}/setup.py build_sphinx