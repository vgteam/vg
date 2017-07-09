#!/usr/bin/env julia

using ArgParse
using Word2Vec

s = ArgParseSettings()
@add_arg_table s begin
    "--model", "-m"
    arg_type = String
    help = "use this model file (write given -i, query given -q)"
    required = true
    default = ""
    "--input", "-i"
    arg_type = String
    help = "input text corpus to model"
    default = ""
    "--kmer-size", "-k"
    arg_type = Int
    help = "kmer size to break input DNA sequences into"
    default = 5
    "--kmer-stride", "-s"
    arg_type = Int
    help = "bases between each successive kmer"
    default = 1
    "--window-size", "-w"
    arg_type = Int
    help = "window size to use during word2vec"
    default = 30
    "--cbow", "-c"
    help = "use continuous bag of words in word2vec"
    action = :store_true
    "--dims", "-d"
    arg_type = Int
    help = "number of dimensions in embedding"
    default = 100
    "--query", "-q"
    arg_type = String
    help = "print the vector for this string"
    default = ""
end

parsed_args = parse_args(ARGS, s)
input = parsed_args["input"]
model_file = parsed_args["model"]
model_size = parsed_args["dims"]
query = parsed_args["query"]
kmer_size = parsed_args["kmer-size"]
kmer_stride = parsed_args["kmer-stride"]
window_size = parsed_args["window-size"]
cbow = parsed_args["cbow"]

function kmers_of(k::Int, j::Int, s::String)
    [s[i:(i+k-1)] for i in range(1,j,Int(floor((length(s)-k)/j))+1)]
end

function kmers_of_file(k::Int, j::Int, f::String)
    kmers = []
    for s in readlines(open(f))
        append!(kmers, kmers_of(k, j, s))
    end
    kmers
end

function write_kmers_of(k::Int, j::Int, input::String, output::String)
    write(open(output, "w"), join(kmers_of_file(k, j, input), " "))
end

if input != ""
    if model_file == ""
        println("an output file is required when building a model")
        exit(1)
    end
    kmers = "$input.kmers"
    write_kmers_of(kmer_size, kmer_stride, input, kmers)
    word2vec(kmers, model_file, verbose = true, size = model_size, window=window_size, cbow=(cbow?1:0))
    rm(kmers)
elseif model_file != ""
    model = wordvectors(model_file)
    if query != ""
        kmers = kmers_of(kmer_size, kmer_stride, query)
        mean_vec = reduce(+, [(haskey(model.vocab_hash, kmer) ? get_vector(model, kmer) : 0) for kmer in kmers])/length(kmers)
        println("$query ", mean_vec)
    end
end
