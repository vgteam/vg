#!/bin/bash
((
jq " .path[]|{(\"\(.name)\"):{\"http://www.w3.org/1999/02/22-rdf-syntax-ns#type\":[{\"value\":\"Path\",\"type\":\"uri\"}]}} " $1|sed '$d';
jq ".edge[]|{(\"\(.from | tostring)\"):{\"before\":[ {\"value\":(\"\(.to |tostring)\"), \"type\" : \"uri\"}]}}" $1 | sed -e "s_^\}_,_"| sed -e 's_^[\{]__'|sed '1d;$d';
jq ".node[]|{(\"\(.id | tostring)\"): {\"http://www.w3.org/1999/02/22-rdf-syntax-ns#value\":[ {\"value\":.sequence, \"type\" : \"literal\"} ] }}" $1|sed -e "s_^\}_,_"| sed -e 's_^[\{]__'|sed '1d;$d';
echo ',';
jq  "  .path[] as \$p|\$p.mapping 
| keys[] as \$k 
| { (\$k | tostring):{\"node\":[ {\"value\":(.[\$k].position.node_id | tostring) 
, \"type\":\"uri\"} ]
, \"step\" : 
	[{\"value\":(\$k | tostring)
	, \"type\":\"literal\"}]
,\"path\":
	[{\"value\":(\$p.name),
	\"type\":\"uri\" } ]}} 
" $1 |sed -e "s_^\}_,_"| sed -e 's_^[\{]__';)|sed '$d';echo '}')
