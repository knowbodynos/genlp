#!/bin/bash

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

sphinx-apidoc -f -o ${DIR}/docs/source/ ${DIR}
python ${DIR}/setup.py build_sphinx