"%PYTHON%" setup.py install
"%PYTHON%" setup.py test
if errorlevel 2 exit 1
