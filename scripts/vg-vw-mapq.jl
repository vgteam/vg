#!/usr/bin/env julia

using ArgParse
using JSON
using Word2Vec

s = ArgParseSettings()
@add_arg_table s begin
    "--vwify", "-v"
    help = "convert JSON-GAM input to vw compatible format"
    action = :store_true
    "--importance", "-i"
    help = "give errors a weight proportional to their mapping quality"
    action = :store_true
    "--phred", "-p"
    help = "regressor already outputs phred format, so don't convert"
    action = :store_true
    "--min-mq", "-q"
    help = "take the minimum of the previous and new mapping quality"
    action = :store_true
    "--mq-model", "-m"
    arg_type = String
    help = "set the MQ using the given model"
    default = ""
    "--dna-embedding", "-e"
    arg_type = String
    help = "use the DNA sequence embedding to generate sequence features"
    default = ""
    "--vw-opts", "-o"
    arg_type = String
    help = "additional options to pass to vw"
    default = ""
    "--kmer-size", "-k"
    arg_type = Int
    help = "kmer size to break input DNA sequences into"
    default = 5
    "--kmer-stride", "-s"
    arg_type = Int
    help = "bases between each successive kmer"
    default = 1
    "--cigar", "-c"
    help = "write the cigar base-by-base in a namespace"
    action = :store_true
    "--sequence", "-S"
    help = "write the read sequence base-by-base"
    action = :store_true
end

parsed_args = parse_args(ARGS, s)

embedding_file = parsed_args["dna-embedding"]
embedding = (embedding_file != "" ? wordvectors(embedding_file) : nothing)
kmer_size = parsed_args["kmer-size"]
kmer_stride = parsed_args["kmer-stride"]
write_cigar = parsed_args["cigar"]
write_sequence = parsed_args["sequence"]

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

function applymodel(alns, model, params)
    (path,io) = mktemp()
    for aln in alns
        println(io, vwify(aln, kmer_size, kmer_stride))
    end
    close(io)
    lines = split(chomp(readstring(pipeline(path, `vw -t -p /dev/stdout -i $model --quiet $(split(params, " "))`))), "\n")
    rm(path)
    assert(length(lines) == length(alns))
    for (line,aln) in zip(lines, alns)
        (result,name) = split(line, " ")
        assert(name == aln["name"])
        result = parsed_args["phred"] ? Int(round(float(result))) : prob2mq(float(result))
        mapping_quality = haskey(aln, "mapping_quality") ? aln["mapping_quality"] : 0
        aln["mapping_quality"] = (parsed_args["min-mq"] ? min(mapping_quality, result) : result)
    end
    return alns
end

function kmers_of(k::Int, j::Int, s::String)
    [s[i:(i+k-1)] for i in range(1,j,Int(floor((length(s)-k)/j))+1)]
end

function seq_embedding(s::String, k::Int, j::Int)
    function in_vocab(kmer::String)
        haskey(embedding.vocab_hash, kmer)
    end
    kmers = filter(in_vocab, kmers_of(k, j, s))
    reduce(+, [get_vector(embedding, kmer) for kmer in kmers])/length(kmers)
end

function vwify(aln, kmer_size::Int, kmer_stride::Int, mqimportance=false)
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
    if write_cigar
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
    end
    sequence = write_sequence ? string("|bases ", join([string(x, y) for (x,y) in zip(1:length(aln["sequence"]),aln["sequence"])], " ")) : ""
    seqvec = (embedding == nothing ? "" : string("|dnavec ", join(["d$i:$(Float16(j))" for (i,j) in enumerate(seq_embedding(aln["sequence"], kmer_size, kmer_stride))], " ")))
    cigar = write_cigar ? string("|cigar ", join(cigar, " ")) : ""
    secondaries = haskey(aln, "secondary_score") ? aln["secondary_score"] : []
    secondaries = join(["s$i:$score" for (i,score) in zip(1:length(secondaries),secondaries)], " ")
    importance = (mqimportance && correct == 0) ? "$(mapqual+2) " : ((mqimportance && correct == 1 && mapqual < 60) ? "$(60-mapqual+1) " : " ")
    string("$correct $importance '$name |alignment identity:$identity score:$score mapqual:$mapqual |scores s0:$score $secondaries $seqvec $sequence $cigar")
end

alns = []
for line in eachline(STDIN)
    aln = JSON.parse(line)
    if parsed_args["vwify"]
        println(vwify(aln, kmer_size, kmer_stride, parsed_args["importance"]))
    elseif parsed_args["mq-model"] != ""
        # we should run the model on the input and set a new mapping quality
        push!(alns, aln)
        if length(alns) >= 1000
            assert(length(alns) > 0)
            for aln in applymodel(alns, parsed_args["mq-model"], parsed_args["vw-opts"])
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
    for aln in applymodel(alns, parsed_args["mq-model"], parsed_args["vw-opts"])
        println(JSON.json(aln))
    end
    alns = []
end
