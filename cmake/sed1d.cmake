file(READ "/Users/zt-home/vg/doc/man/vg-manpage.md" _content)
string(REGEX REPLACE "^[^\n]*\n" "" _stripped "${_content}")
file(WRITE "/Users/zt-home/vg/doc/wiki/vg-manpage.md" "${_stripped}")
