# We used to have a script here to set up all the include and library search
# paths for the vg build. But now the Makefile knows how to do it all for the
# build, and the vg binary knows where to look for its dynamic libraries.
echo 1>&2 "Sourcing source_me.sh is no longer necessary"
