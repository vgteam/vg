.PHONY: clean pre

vg.o: $(PROTO_OUT_DIR)/vg.pb.h $(LIB_VCFLIB) $(FASTAHACK_DIR)/Fasta.o $(LIB_GSSW) $(LIB_SPARSEHASH) $(LRU_CACHE_DIR)/lru_cache.h stream.hpp $(LIB_PROTOBUF) $(LIB_SDSLLITE)
        $(CXX) $(CXXFLAGS) -c -o vg.o vg.cpp $(INCLUDES)

gssw_aligner.o: gssw_aligner.cpp gssw_aligner.hpp cpp/vg.pb.h $(LIB_GSSW) $(LIB_PROTOBUF) $(LIB_SPARSEHASH)
        $(CXX) $(CXXFLAGS) -c -o gssw_aligner.o gssw_aligner.cpp $(INCLUDES)

vg_set.o: vg_set.cpp vg_set.hpp vg.hpp index.hpp cpp/vg.pb.h $(LIBGSSW) $(LIBPROTOBUF) $(SPARSEHASH) $(SDSLLITE)
        $(CXX) $(CXXFLAGS) -c -o vg_set.o vg_set.cpp $(INCLUDES)

mapper.o: mapper.cpp mapper.hpp cpp/vg.pb.h $(LIBPROTOBUF) $(SPARSEHASH) $(SDSLLITE)
        $(CXX) $(CXXFLAGS) -c -o mapper.o mapper.cpp $(INCLUDES)

main.o: main.cpp $(LIBVCFLIB) $(fastahack/Fasta.o) $(LIBGSSW) stream.hpp  $(LIBPROTOBUF) $(SPARSEHASH) $(SDSLLITE)
        $(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)


region.o: region.cpp region.hpp $(LIBPROTOBUF) $(SPARSEHASH)
        $(CXX) $(CXXFLAGS) -c -o region.o region.cpp $(INCLUDES)

index.o: index.cpp index.hpp $(LIBPROTOBUF) $(SPARSEHASH)
        $(CXX) $(CXXFLAGS) -c -o index.o index.cpp $(INCLUDES)

utility.o: utility.cpp utility.hpp $(LIBPROTOBUF) $(SPARSEHASH)
        $(CXX) $(CXXFLAGS) -c -o utility.o utility.cpp $(INCLUDES)

path.o: path.cpp path.hpp $(LIBPROTOBUF) $(SPARSEHASH)
        $(CXX) $(CXXFLAGS) -c -o path.o path.cpp $(INCLUDES)

edit.o: edit.cpp edit.hpp $(LIBPROTOBUF)
        $(CXX) $(CXXFLAGS) -c -o edit.o edit.cpp $(INCLUDES)

alignment.o: alignment.cpp alignment.hpp $(LIBHTS)  $(LIBPROTOBUF) $(SPARSEHASH) edit.hpp edit.cpp
        $(CXX) $(CXXFLAGS) -c -o alignment.o alignment.cpp $(INCLUDES)

sha1/sha1.o: sha1/sha1.cpp sha1/sha1.hpp
        $(CXX) $(CXXFLAGS) -c -o sha1/sha1.o sha1/sha1.cpp $(INCLUDES)

json2pb.o: json2pb.cpp json2pb.h bin2ascii.h $(LIBPROTOBUF)
        $(CXX) $(CXXFLAGS) -c -o json2pb.o json2pb.cpp $(INCLUDES)

entropy.o: entropy.cpp entropy.hpp
        $(CXX) $(CXXFLAGS) -c -o entropy.o entropy.cpp $(INCLUDES)

pileup.o: pileup.cpp pileup.hpp cpp/vg.pb.h vg.hpp stream.hpp json2pb.h $(LIBPROTOBUF) $(SPARSEHASH)
        $(CXX) $(CXXFLAGS) -c -o pileup.o pileup.cpp $(INCLUDES)

caller.o: caller.cpp caller.hpp cpp/vg.pb.h vg.hpp stream.hpp json2pb.h pileup.hpp $(LIBPROTOBUF) $(SPARSEHASH)
        $(CXX) $(CXXFLAGS) -c -o caller.o caller.cpp $(INCLUDES)


