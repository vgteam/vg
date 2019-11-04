// #include <iostream>

// #include <seqan/basic.h>
// #include <seqan/stream.h>   // For printing strings.
// #include <seqan/score.h>  // The module score.

// using namespace seqan;
// // Extend SeqAn by a user-define scoring matrix.
// namespace seqan {

// // We have to create a new specialization of the ScoringMatrix_ class
// // for the DNA alphabet.  For this, we first create a new tag.
// struct UserDefinedMatrix {};

// // Then, we specialize the class ScoringMatrix_ for the Dna5 alphabet.
// template <>
// struct ScoringMatrixData_<int, Dna5, UserDefinedMatrix>
// {
//     enum
//     {
//         VALUE_SIZE = ValueSize<Dna5>::VALUE,
//         TAB_SIZE = VALUE_SIZE * VALUE_SIZE
//     };

//     static inline int const * getData()
//     {
//         // The user defined data table.  In this case, we use the data from BLOSUM-30.
//         static int const _data[TAB_SIZE] =
//         {
//             1, 0, 0, 0, 0,
//             0, 1, 0, 0, 0,
//             0, 0, 1, 0, 0,
//             0, 0, 0, 1, 0,
//             0, 0, 0, 0, 0
//         };
//         return _data;
//     }

// };
// }  // namespace seqan
// We define a function showScoringMatrix for displaying a matrix.

// // Print a scoring scheme matrix to stdout.
// template <typename TScoreValue, typename TSequenceValue, typename TSpec>
// void showScoringMatrix(Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & scoringScheme)
// {
//     // Print top row.
//     for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
//         std::cout << "\t" << TSequenceValue(i);
//     std::cout << std::endl;
//     // Print each row.
//     for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
//     {
//         std::cout << TSequenceValue(i);
//         for (unsigned j = 0; j < ValueSize<TSequenceValue>::VALUE; ++j)
//         {
//             std::cout << "\t" << score(scoringScheme, TSequenceValue(i), TSequenceValue(j));
//         }
//         std::cout << std::endl;
//     }
// }
// int main()
// {
//     // 1. Define type and constants.
//     //
//     // Define types for the score value and the scoring scheme.
//     typedef int TValue;
//     typedef Score<TValue, ScoreMatrix<Dna5, Default> > TScoringScheme;
//     // Define our gap scores in some constants.
//     int const gapOpenScore = -1;
//     int const gapExtendScore = -1;

//     // 2. Construct scoring scheme with default/empty matrix.
//     //
//     // Construct new scoring scheme, alternatively only give one score
//     // that is used for both opening and extension.
//     TScoringScheme scoringScheme(gapExtendScore, gapOpenScore);

//     // 3. Fill the now-existing ScoreMatrix
//     //
//     // The scoring scheme now already has a matrix of the size
//     // ValueSize<Dna5>::VALUE x ValueSize<Dna5>::VALUE which
//     // we can now fill.

//     // 3.1 We fill the scoring scheme with the product of the coordinates.
//     std::cout << std::endl << "Coordinate Products" << std::endl;
//     for (unsigned i = 0; i < ValueSize<Dna5>::VALUE; ++i)
//     {
//         for (unsigned j = 0; j < ValueSize<Dna5>::VALUE; ++j)
//         {
//             setScore(scoringScheme, Dna5(i), Dna5(j), i * j);
//         }
//     }
//     showScoringMatrix(scoringScheme);

//     // 3.2 Now, we fill it with the user defined matrix above.
//     std::cout << "User defined matrix (also Dna5 scoring matrix)..." << std::endl;
//     setDefaultScoreMatrix(scoringScheme, UserDefinedMatrix());
//     showScoringMatrix(scoringScheme);

//     // 4. Show our user-defined Dna5 scoring matrix.
//     std::cout << "User DNA scoring scheme..." << std::endl;
//     Score<TValue, ScoreMatrix<Dna5, UserDefinedMatrix> > userScoringSchemeDna;
//     showScoringMatrix(userScoringSchemeDna);

//     return 0;
// }
// Here is the output of the program:

// Coordinate Products
// 	A	C	G	T	N
// A	0	0	0	0	0
// C	0	1	2	3	4
// G	0	2	4	6	8
// T	0	3	6	9	12
// N	0	4	8	12	16
// User defined matrix (also Dna5 scoring matrix)...
// 	A	C	G	T	N
// A	1	0	0	0	0
// C	0	1	0	0	0
// G	0	0	1	0	0
// T	0	0	0	1	0
// N	0	0	0	0	0
// User DNA scoring scheme...
// 	A	C	G	T	N
// A	1	0	0	0	0
// C	0	1	0	0	0
// G	0	0	1	0	0
// T	0	0	0	1	0
// N	0	0	0	0	0
// Loading Score Matrices From File
// This small demo program shows how to load a score matrix from a file. Examples for score file are demos/howto/scores/dna_example.txt for DNA alphabets and tests/sPAM250 for amino acids.

// Include the necessary headers.

// #include <iostream>

// #include <seqan/basic.h>
// #include <seqan/stream.h>   // For printing strings.
// #include <seqan/score.h>  // The module score.

// using namespace seqan;
// We define a function that can show a scoring matrix.

// // Print a scoring scheme matrix to stdout.
// template <typename TScoreValue, typename TSequenceValue, typename TSpec>
// void showScoringMatrix(Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & scoringScheme)
// {
//     // Print top row.
//     for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
//         std::cout << "\t" << TSequenceValue(i);
//     std::cout << std::endl;
//     // Print each row.
//     for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
//     {
//         std::cout << TSequenceValue(i);
//         for (unsigned j = 0; j < ValueSize<TSequenceValue>::VALUE; ++j)
//         {
//             std::cout << "\t" << score(scoringScheme, TSequenceValue(i), TSequenceValue(j));
//         }
//         std::cout << std::endl;
//     }
// }
// Finally, the main program loads the scoring matrix and then shows it.

// int main(int argc, char ** argv)
// {
//     typedef int TScoreValue;

//     Score<TScoreValue, ScoreMatrix<Dna5, Default> > scoreMatrix;
//     loadScoreMatrix(scoreMatrix, toCString(getAbsolutePath("demos/howto/scores/dna_example.txt")));
//     showScoringMatrix(scoreMatrix);

//     return 0;
// }
// Here’s the program output.

//     A   C   G   T
// A   1   -1  -1  -1
// C   -1  1   -1  -1
// G   -1  -1  1   -1
// T   -1  -1  -1  1
// © Copyright 2015, The SeqAn Team. Revision 88e3a0bb.

// Built with Sphinx using a theme provided by Read the Docs.