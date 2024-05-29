cp -f conf.py src_build/
cd src_build
sphinx-build -a -E . ../_build
