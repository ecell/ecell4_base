set "CMAKE_GENERATOR=Visual Studio 15 2017"
"%PYTHON%" setup.py install
"%PYTHON%" setup.py test
if errorlevel 1 exit 1
