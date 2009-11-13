#!/bin/sh
# Run this to generate all the initial makefiles, etc.
# This was lifted from the Gimp, and adapted slightly by
# Raph Levien .

DIE=0

PROJECT="brown"

(autoconf --version) < /dev/null > /dev/null 2>&1 || {
        echo
        echo "You must have autoconf installed to compile $PROJECT."
        DIE=1
}

(libtoolize --version) < /dev/null > /dev/null 2>&1 || {
        echo
        echo "You must have libtool installed to compile $PROJECT."
        DIE=1
}

(automake --version) < /dev/null > /dev/null 2>&1 || {
        echo
        echo "You must have automake installed to compile $PROJECT."
        DIE=1
}

if test "$DIE" -eq 1; then
        exit 1
fi

libtoolize -c --force --automake

case $CC in
*xlc | *xlc\ * | *lcc | *lcc\ *) am_opt=--include-deps;;
esac

for dir in .
  do 
  echo -n Running autotools for $dir ...
  (cd $dir; \
  { echo -n ' aclocal '; aclocal; } && \
  { echo -n 'autoheader '; autoheader -f ; } && \
  { echo -n 'automake ';  automake --copy --add-missing $am_opt; } && \
  { echo -n 'autoconf '; autoconf; } && \
  echo )
  
  if  test $? != 0 ; then
      echo "Error processing $dir"
      exit $? 
  fi

done

echo 'Finished running autotools.  Run ./configure next.'


