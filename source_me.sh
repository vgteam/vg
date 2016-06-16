export LIBRARY_PATH=`pwd`/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
export LD_INCLUDE_PATH=`pwd`/include:$LD_INCLUDE_PATH
export C_INCLUDE_PATH=`pwd`/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=`pwd`/include:$CPLUS_INCLUDE_PATH
export INCLUDE_PATH=`pwd`/include:$INCLUDE_PATH
export PATH=`pwd`/bin:`pwd`/scripts:$PATH
export CC=$(which gcc)
export CXX=$(which g++)

# add bash autocompletion
if [ -n "$BASH" ]
then
	 _vg_complete()
	 {
		  local cur prev
		  COMPREPLY=()
		  cur="${COMP_WORDS[COMP_CWORD]}"
		  opts="construct deconstruct view vectorize index find paths align map stats join ids concat kmers sim mod surject msga pileup call genotype compare scrub circularize validate version"

		  COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
		  return 0
	 }
	 complete -o nospace -F _vg_complete vg
fi