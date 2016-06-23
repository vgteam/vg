_vg_complete()
{
	 local cur prev
	 COMPREPLY=()
	 cur="${COMP_WORDS[COMP_CWORD]}"
	 opts="construct deconstruct view vectorize index find paths align map stats join ids concat kmers sim mod surject msga pileup call genotype compare scrub circularize validate version"
	 
	 if [ $COMP_CWORD -eq 1 ]
	 then
		  COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
	 else
		  COMPREPLY=( $(compgen -W "$(ls)" -- $cur) )
	 fi
	 return 0
}
complete -o nospace -F _vg_complete vg
