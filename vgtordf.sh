#!/bin/bash
((
jq " .path[]|{(\"\(.name)\"):{\"rdf:type\":[{\"value\":\"Path\",\"type\":\"uri\"}]}} " $1|sed '$d';
jq ".edge[]|{(\"\(.from | tostring)\"):{\"before\":[ {\"value\":(\"\(.to |tostring)\"), \"type\" : \"uri\"}]}}" $1 | sed -e "s_^\}_,_"| sed -e 's_^[\{]__'|sed '1d;$d';
jq ".node[]|{(\"\(.id | tostring)\"): {\"rdf:value\":[ {\"value\":.sequence, \"type\" : \"literal\"} ] }}" $1|sed -e "s_^\}_,_"| sed -e 's_^[\{]__'|sed '1d;$d';
echo ',';
jq  "  .path[].mapping | keys[] as \$k | { (\$k | tostring):{\"node\":[ {\"value\":(.[\$k].position.node_id | tostring) , \"type\":\"uri\"} ], \"step\" : [{\"value\":(\$k | tostring), \"type\":\"literal\"} ]} }" $1 |sed -e "s_^\}_,_"| sed -e 's_^[\{]__';)|sed '$d';echo '}')
