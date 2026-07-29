_vg_complete()
{
	 local cur prev
	 COMPREPLY=()
	 cur="${COMP_WORDS[COMP_CWORD]}"
	 opts="add align annotate augment autoindex bench-dist-query benchmark call chain chains chunk circularize clip cluster combine concat construct convert deconstruct depth describe dotplot filter find gamcompare gampcompare gamsort gbwt genotype giraffe haplotypes ids index inject kmers map mask mcmc minimizer mod mpmap pack paths primers prune rna sim simplify snarls sort stats surject trace translate validate vectorize version view viz zipcode"
	 
	 if [ $COMP_CWORD -eq 1 ]
	 then
		  COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
	 else
		  COMPREPLY=( $(compgen -W "$(ls)" -- $cur) )
	 fi
	 return 0
}
complete -o nospace -F _vg_complete vg
