#!/usr/bin/env julia

using ArgParse
using JSON

s = ArgParseSettings()
@add_arg_table s begin
    "--vwify", "-v"
    help = "convert JSON-GAM input to vw compatible format"
    action = :store_true
    "--mq-model", "-m"
    arg_type = String
    help = "set the MQ using the given model"
    default = ""
end

parsed_args = parse_args(ARGS, s)

global MQ_MAX = 60

function float2phred(prob)
    return min(-10 * log10(prob), MQ_MAX)
end

function phred2float(qual)
    return 10 ^ (qual * -.1)
end

function prob2mq(prob)
    v = 1-prob
    if v == 1
        return 0
    elseif v > 0
        return min(MQ_MAX, Int32(round(float2phred(v))))
    else
        return MQ_MAX
    end
end

function applymodel(alns, model)
    (path,io) = mktemp()
    for aln in alns
        println(io, vwify(aln))
    end
    close(io)
    lines = split(chomp(readstring(pipeline(path, `vw -t -p /dev/stdout -i $model --quiet`))), "\n")
    rm(path)
    assert(length(lines) == length(alns))
    for (line,aln) in zip(lines, alns)
        (result,name) = split(line, " ")
        assert(name == aln["name"])
        aln["mapping_quality"] = prob2mq(float(result))
    end
    return alns
end

function vwify(aln)
    name = haskey(aln, "name") ? aln["name"] : ""
    correct = haskey(aln, "correct") ? aln["correct"] : 0
    correct = correct > 0.5 ? 1 : 0
    if haskey(aln, "correct") aln["correct"] = correct end
    score = haskey(aln, "score") ? aln["score"] : 0
    identity = haskey(aln, "identity") ? aln["identity"] : 0
    mapqual = haskey(aln, "mapping_quality") ? aln["mapping_quality"] : 0
    path = haskey(aln, "path") && haskey(aln["path"], "mapping") ? aln["path"]["mapping"] : []
    i = 1
    cigar = []
    for mapping in path
        for elem in mapping["edit"]
            #println(elem)
            fromlen = haskey(elem, "from_length") ? elem["from_length"] : 0
            tolen = haskey(elem, "to_length") ? elem["to_length"] : 0
            seq = haskey(elem, "sequence") ? elem["sequence"] : ""
            if fromlen == tolen
                if seq != ""
                    for j in range(i, tolen) push!(cigar, string(j, "S")) end
                else
                    for j in range(i, tolen) push!(cigar, string(j, "M")) end
                end
            elseif fromlen > tolen
                for j in range(i, fromlen) push!(cigar, string(j, "D")) end
            elseif tolen > fromlen
                for j in range(i, tolen) push!(cigar, string(j, "I")) end
            end
            i += tolen
        end
    end
    sequence = join([string(x, y) for (x,y) in zip(1:length(aln["sequence"]),aln["sequence"])], " ")
    cigar = join(cigar, " ")
    secondaries = haskey(aln, "secondary_score") ? aln["secondary_score"] : []
    secondaries = join(["s$i:$score" for (i,score) in zip(1:length(secondaries),secondaries)], " ")
    string(correct, " '", name, " |alignment identity:$identity score:$score mapqual:$mapqual", " |others ", secondaries, " |sequence ", sequence, " |cigar ", cigar)
end

alns = []
for line in eachline(STDIN)
    aln = JSON.parse(line)
    if parsed_args["vwify"]
        println(vwify(aln))
    elseif parsed_args["mq-model"] != ""
        # we should run the model on the input and set a new mapping quality
        push!(alns, aln)
        if length(alns) >= 1000
            assert(length(alns) > 0)
            for aln in applymodel(alns, parsed_args["mq-model"])
                println(JSON.json(aln))
            end
            alns = []
        end
    else
        # default is to do nothing to input
        println(JSON.json(aln))
    end
end

# clean up the buffer
if parsed_args["mq-model"] != "" && length(alns) > 0
    for aln in applymodel(alns, parsed_args["mq-model"])
        println(JSON.json(aln))
    end
    alns = []
end
