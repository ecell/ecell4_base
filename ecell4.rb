require 'formula'

# Contact celery [at] g-language.org for questions

# Or should I call this formula plain ecell...?
class Ecell4 < Formula
  homepage 'http://www.e-cell.org/'
  url 'https://github.com/ecell/ecell4/tarball/master'
  sha1 '2b1bdf70b0bdfd09e091943eb110b1aab03db6c6'
  version '4'

  # Bad practice for finding path to "cython"
  tmp = `mdfind -name cython | grep /bin/ | head -n 1`.split("/")
  tmp.pop
  CYTHONPATH = tmp.join("/")

  # Define some constants for environment variables
  PATH = String.new(ENV['PATH']) + ":" + CYTHONPATH
  PREFIX = "/usr/local/Cellar/ecell4/4"
  tmp = /\d\.\d\.\d/.match(`python --version 2>&1`).to_s.split(".")
  tmp.pop
  PYTHONVERSION = tmp.join(".")
  PYTHONPATH = PREFIX + "/lib/python" + PYTHONVERSION + "/site-packages"
  NEWPATH = "/usr/local/lib/python" + PYTHONVERSION + "/site-packages"
  CPATH = PREFIX + "/include"
  LIBRARY_PATH = PREFIX + "/lib"
  tmp = `echo $SHELL`.split("/")
  SHELL = tmp.pop.chomp

  # Uncomment lines to include in installation
  # NOTE - Some of these modules may require additional Python modules
  targets = ["core",       "core_python",
             # "egfrd",      "egfrd_python", # Requires SciPy
             "gillespie",  "gillespie_python",
             "ode",        "ode_python",
             # "spatiocyte", "spatiocyte_python", # Requires ecs
             "reaction_reader"]

  # Dependencies
  depends_on 'pkg-config'
  depends_on 'gsl'
  depends_on 'boost'
  depends_on 'hdf5' => '--enable-cxx'

  # This will automatically build the target module
  def waf(target)
    puts "Building " + target + "..."

    tmp = `mdfind -name cython | grep /bin/ | head -n 1`.split("/")
    tmp.pop
    tmp = tmp.join("/")

    ENV['PATH'] = PATH + ":" + tmp
    ENV['CPATH'] = CPATH
    ENV['LIBRARY_PATH'] = LIBRARY_PATH
    ENV['PYTHONPATH'] = PYTHONPATH

    # Process the waf stuff
    Dir.chdir(target)
    system "../waf configure --prefix=" + PREFIX
    system "../waf build"
    system "../waf install"
    Dir.chdir("../")
  end

  def install
    ENV['PATH'] = PATH

    # Install pip if not present
    print "Looking for pip..."
    unless `which pip`.length > 1
      puts " not found"
      puts "\n\e[31mWarning:\e[0m"
      puts "We are trying to install \e[32mpip\e[0m into your environment"
      puts "If you do not wish for homebrew to perform this task, exit"
      puts "and install pip by performing"
      puts "\n  $ sudo easy_install pip"
      puts "\nIf you wish to proceed, enter root password"
      system "sudo easy_install pip"
    else
      puts " found"
    end

    # Install cython if not present
    print "Looking for cython..."
    if `pip freeze | grep -ic ^cython`.chomp.to_i == 0
      puts " not found"
      puts "\n\e[31mWarning:\e[0m"
      puts "We are trying to install \e[32mcython\e[0m into your environment"
      puts "If you do not wish for homebrew to perform this task, exit" 
      puts "and install cython by performing"
      puts "\n  $ sudo pip install cython"
      puts "\nor if you are not root"
      puts "\n  $ pip install cython --user"
      puts "\nIf you wish to proceed, enter root password"
      system "sudo pip install cython"
    else
      puts " found"
    end

    # Install scipy if not present
    print "Looking for scipy..."
    if `pip freeze | grep -ic ^scipy`.chomp.to_i == 0
      puts " not found: egfrd will not be installed"
      # puts " not found"
      # system "sudo pip install scipy"
    else
      puts " skipped: egfrd will not be installed" # Build error
      # puts " found"
      # targets.push("egfrd", "egfrd_python")
    end

    # Install ecs if not present
    print "Looking for ecs..."
    if `pip freeze | grep -ic ^ecs`.chomp.to_i == 0
      puts " not found: spatiocyte will not be installed"
      # puts " not found"
      # system "sudo pip install ecs"
    else
      puts " found: installing spatiocyte"
      targets.push("spatiocyte", "spatiocyte_python")
    end
    
    targets.each {
      |target|
      waf(target)
    }

    if SHELL == "tcsh" || SHELL == "csh"
      puts "Remember to add"
      puts "'setenv PYTHONPATH $PYTHONPATH:#{NEWPATH}'"
      puts "to your .#{SHELL}rc"
    else
      puts "Remember to add"
      puts "'export PYTHONPATH=$PYTHONPATH:#{NEWPATH}'"
      puts "to your .#{SHELL}rc"
    end
  end

  test do
    if SHELL == "tcsh" || SHELL == "csh"
      system "setenv PYTHONPATH $PYTHONPATH:#{NEWPATH}"
    else
      system "export PYTHONPATH=$PYTHONPATH:#{NEWPATH}"
    end

    system "python -c 'import ecell4'"
  end
end
