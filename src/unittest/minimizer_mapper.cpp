/// \file minimizer_mapper.cpp
///  
/// unit tests for the minimizer mapper

#include <iostream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../minimizer_mapper.hpp"
#include "../build_index.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

// We define a child class to expose protected stuff for testing
class TestMinimizerMapper : public MinimizerMapper {
public:
    TestMinimizerMapper(
        gbwtgraph::GBWTGraph gbwt_graph,
        vector<gbwtgraph::DefaultMinimizerIndex*> minimizer_indexes,
        MinimumDistanceIndex distance_index,
        PathPositionHandleGraph* handle_graph) 
        : MinimizerMapper(gbwt_graph, minimizer_indexes, distance_index, handle_graph){};
    using MinimizerMapper::MinimizerMapper;
    using MinimizerMapper::score_extension_group;
    using MinimizerMapper::Minimizer;
    using MinimizerMapper::fragment_length_distr;
};

TEST_CASE("MinimizerMapper::score_extension_group works", "[giraffe][mapping]") {

    // Define an Alignment to score extensions against.
    Alignment aln;
    aln.set_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    
    // Define a bunch of fake GaplessExtensions to chain up and score
    vector<GaplessExtension> to_score;
    
    SECTION("Score of no gapless extensions is 0") {
        REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 0);
    }
    
    SECTION("One-base extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 1;
        to_score.back().read_interval.second = 2;
        to_score.back().score = 1;
        
        SECTION("Score of a 1-base gapless extension is 1") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 1);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 2;
        to_score.back().read_interval.second = 3;
        to_score.back().score = 1;
        
        SECTION("Score of two adjacent 1-base gapless extensions is 2") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 2);
        }
        
    }
    
    SECTION("Longer extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 1;
        to_score.back().read_interval.second = 10;
        to_score.back().score = 9;
        
        SECTION("Score of one 9-base extension is 9") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 9);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 10;
        to_score.back().read_interval.second = 19;
        to_score.back().score = 9;
        
        SECTION("Score of two 9-base extensions abutting is 18") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 18);
        }
        
        to_score.back().read_interval.first++;
        to_score.back().read_interval.second++;
        
        SECTION("Score of two 9-base extensions separated by a 1-base gap is 12") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 12);
        }
        
        to_score.back().read_interval.first++;
        to_score.back().read_interval.second++;
        
        SECTION("Score of two 9-base extensions separated by a 2-base gap is 11") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 11);
        }
        
        to_score.back().read_interval.first -= 3;
        to_score.back().read_interval.second -= 3;
        
        SECTION("Score of two 9-base extensions with a 1-base ovealap is 12") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 11);
        }
    }
    
    SECTION("Many possibly overlapping extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 0;
        to_score.back().read_interval.second = 1;
        to_score.back().score = 1;
    
        for (size_t i = 0; i < 35; i++) {
        
            to_score.emplace_back();
            to_score.back().read_interval.first = i + 1;
            to_score.back().read_interval.second = i + 1 + 30;
            to_score.back().score = 30;
        
        }
    
        
        
        SECTION("Score of one 1-base extensions and 2 30-base extensions is 61") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 61);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 28;
        to_score.back().read_interval.second = 28 + 45;
        to_score.back().score = 45;
        
        // Sort by read interval as is required
        std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
            return (a.read_interval.first < b.read_interval.first) ||
                (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
        });
        
        SECTION("Score of one 1-base extension, a 30-base extension, a backtrack of 4, and a 45-base extension is correct") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == (31 + 45 - 6 - 3));
        }
        
        
        for (size_t i = 3; i < 29; i++) {
             to_score.emplace_back();
            to_score.back().read_interval.first = i;
            to_score.back().read_interval.second = i + 15;
            to_score.back().score = 15;
        }
        
        // Sort by read interval as is required
        std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
            return (a.read_interval.first < b.read_interval.first) ||
                (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
        });
        
        SECTION("Score of one 1-base extension, a 30-base extension, a backtrack of 4, and a 45-base extension is not distracted") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == (31 + 45 - 6 - 3));
        }
    
    }
}

TEST_CASE("Fragment length distribution gets reasonable value", "[giraffe][mapping]") {

        vector<int64_t> distances { 27, 69, 76, 88, 107, 114, 119, 121, 124, 124, 125, 125, 125, 125, 126, 126, 126, 127, 127, 127, 127, 127, 127, 
            128, 128, 128, 128, 128, 129, 129, 129, 129, 130, 130, 131, 131, 131, 131, 132, 132, 132, 133, 133, 133, 133, 133, 134, 134, 134, 134,
            134, 134, 134, 134, 135, 135, 136, 136, 136, 136, 136, 136, 136, 136, 137, 137, 137, 138, 139, 139, 139, 140, 140, 140, 140, 140, 140, 
            141, 141, 141, 141, 141, 141, 141, 141, 141, 142, 142, 142, 142, 142, 142, 143, 143, 143, 143, 143, 143, 143, 144, 144, 144, 144, 144, 
            144, 144, 144, 145, 145, 145, 145, 145, 145, 145, 145, 145, 146, 146, 146, 146, 146, 147, 147, 147, 147, 147, 148, 148, 148, 148, 148, 
            148, 149, 149, 149, 149, 150, 150, 150, 151, 151, 151, 151, 151, 151, 152, 152, 152, 153, 153, 153, 153, 153, 153, 153, 153, 154, 154, 
            154, 154, 154, 154, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 156, 156, 156, 156, 156, 156, 156, 156, 
            156, 156, 157, 157, 157, 157, 158, 158, 158, 158, 158, 158, 159, 159, 159, 159, 159, 159, 159, 160, 160, 160, 160, 160, 160, 160, 160, 
            160, 160, 160, 160, 160, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 162, 162, 162, 162, 163, 163, 163, 
            164, 164, 164, 164, 164, 164, 164, 165, 165, 165, 165, 165, 165, 165, 165, 166, 166, 166, 166, 166, 166, 166, 166, 166, 166, 167, 167, 
            167, 167, 167, 167, 167, 167, 167, 168, 168, 168, 168, 169, 169, 169, 169, 169, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 171, 
            171, 171, 171, 171, 171, 171, 172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 172, 173, 173, 173, 173, 173, 173, 173, 173, 
            174, 174, 174, 174, 174, 175, 175, 175, 175, 175, 175, 176, 176, 176, 176, 176, 176, 176, 177, 178, 178, 178, 178, 178, 178, 178, 178,
            178, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 180, 180, 180, 180, 180, 180, 180, 181, 181, 181, 181, 181, 181, 181, 181,
            182, 182, 182, 182, 182, 182, 182, 182, 183, 183, 183, 183, 183, 183, 183, 183, 183, 184, 184, 184, 184, 184, 185, 185, 185, 185, 185,
            186, 186, 186, 186, 186, 186, 186, 187, 187, 187, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 189, 189, 189,
            189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 191, 191, 191, 
            192, 192, 192, 192, 192, 192, 192, 193, 193, 193, 193, 193, 193, 193, 194, 194, 194, 194, 194, 194, 194, 194, 195, 195, 195, 195, 195, 
            195, 195, 195, 195, 195, 196, 196, 196, 196, 196, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 198, 198, 198, 198, 198, 198, 198, 
            198, 198, 198, 199, 199, 199, 199, 199, 199, 199, 200, 200, 200, 200, 200, 200, 200, 200, 201, 201, 201, 201, 201, 201, 201, 201, 201, 
            202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 203, 203, 203, 203, 203, 203, 203, 203, 204, 204, 204, 204, 204, 204, 204, 204, 204, 
            204, 204, 204, 204, 205, 205, 205, 205, 205, 205, 205, 205, 205, 206, 206, 206, 206, 206, 206, 206, 206, 207, 207, 207, 207, 207, 207, 
            207, 207, 207, 208, 208, 208, 208, 208, 208, 208, 208, 208, 209, 209, 209, 209, 209, 209, 210, 210, 210, 210, 210, 210, 210, 211, 211, 
            211, 211, 211, 211, 211, 212, 212, 212, 212, 212, 212, 212, 213, 213, 213, 213, 213, 213, 213, 213, 214, 214, 214, 214, 214, 214, 214, 
            214, 214, 214, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 216, 216, 216, 216, 216, 216, 216, 216, 216, 217, 217, 217, 217, 
            217, 217, 217, 217, 218, 218, 218, 218, 218, 218, 218, 218, 218, 218, 218, 219, 219, 219, 219, 219, 219, 219, 219, 219, 219, 220, 220, 
            220, 220, 220, 220, 220, 220, 220, 220, 220, 220, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 222, 222, 222, 222, 222, 
            222, 222, 222, 222, 222, 223, 223, 223, 223, 223, 223, 224, 224, 224, 224, 224, 224, 224, 225, 225, 225, 225, 225, 225, 225, 225, 225, 
            226, 226, 226, 226, 226, 226, 226, 226, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 228, 228, 228, 229, 229, 229, 229, 
            229, 230, 230, 230, 230, 230, 230, 230, 230, 231, 231, 231, 231, 231, 231, 231, 231, 232, 232, 232, 232, 232, 232, 232, 233, 233, 233, 
            233, 233, 233, 233, 233, 233, 233, 233, 233, 234, 234, 234, 234, 234, 234, 234, 234, 234, 234, 234, 235, 235, 235, 235, 235, 236, 236, 
            236, 236, 236, 236, 236, 236, 236, 236, 237, 237, 237, 237, 237, 237, 237, 238, 238, 238, 238, 238, 238, 238, 238, 238, 238, 238, 238, 
            238, 238, 238, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 240, 240, 240, 240, 240, 240, 240, 241, 241, 
            241, 241, 241, 241, 241, 241, 241, 241, 242, 242, 242, 242, 242, 242, 242, 243, 243, 243, 243, 243, 243, 243, 244, 244, 244, 244, 244, 
            244, 244, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 246, 246, 246, 246, 246, 246, 247, 247, 247, 
            247, 247, 247, 247, 247, 247, 247, 247, 247, 248, 248, 248, 248, 248, 248, 248, 248, 248, 248, 248, 248, 248, 249, 249, 249, 249, 249, 
            249, 249, 250, 250, 250, 250, 250, 250, 250, 250, 251, 251, 251, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252, 252, 252, 252, 252, 
            253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 255, 255, 255, 255, 
            255, 255, 255, 255, 255, 255, 255, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 257, 257, 257, 257, 257, 257, 257, 257, 
            257, 257, 258, 258, 258, 258, 259, 259, 259, 259, 259, 259, 260, 260, 260, 260, 260, 260, 260, 260, 260, 260, 260, 260, 261, 261, 261, 
            261, 261, 261, 261, 261, 261, 262, 262, 262, 262, 262, 262, 262, 262, 262, 262, 262, 263, 263, 263, 263, 264, 264, 264, 264, 264, 264, 
            264, 264, 264, 264, 265, 265, 265, 265, 265, 265, 265, 265, 265, 266, 266, 266, 266, 266, 266, 266, 266, 266, 267, 267, 267, 268, 268, 
            268, 268, 268, 268, 268, 268, 269, 269, 269, 269, 269, 269, 269, 269, 269, 269, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 
            271, 271, 271, 271, 271, 271, 271, 271, 271, 271, 272, 272, 272, 272, 272, 273, 273, 273, 273, 273, 273, 273, 274, 274, 274, 274, 274, 
            274, 274, 275, 275, 275, 276, 276, 276, 276, 276, 276, 276, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 278, 278, 278, 278, 
            278, 278, 278, 278, 279, 279, 279, 279, 279, 279, 279, 279, 280, 280, 280, 280, 280, 280, 280, 281, 281, 281, 281, 281, 281, 281, 281, 
            281, 281, 281, 281, 281, 281, 281, 281, 282, 282, 282, 282, 282, 282, 282, 282, 283, 283, 283, 283, 283, 283, 283, 283, 284, 284, 284, 
            284, 284, 285, 285, 285, 285, 285, 285, 285, 285, 285, 285, 286, 286, 286, 286, 286, 286, 287, 287, 287, 287, 287, 287, 287, 287, 288,
            288, 288, 288, 288, 288, 288, 289, 289, 289, 289, 289, 289, 289, 290, 290, 290, 291, 291, 291, 291, 291, 291, 291, 291, 291, 292, 292, 
            292, 292, 292, 292, 292, 292, 293, 293, 293, 293, 293, 293, 293, 293, 294, 294, 294, 294, 294, 294, 294, 294, 295, 295, 295, 295, 295,
            295, 295, 296, 296, 296, 296, 296, 296, 296, 296, 297, 297, 297, 297, 297, 297, 297, 298, 298, 298, 298, 298, 298, 298, 298, 298, 299,
            299, 299, 299, 299, 299, 299, 299, 299, 299, 300, 300, 300, 300, 300, 300, 300, 301, 301, 301, 301, 301, 301, 301, 301, 301, 302, 302,
            302, 302, 302, 302, 302, 302, 303, 303, 303, 303, 303, 303, 303, 303, 304, 304, 304, 304, 304, 304, 304, 305, 305, 305, 305, 306, 306,
            306, 306, 307, 307, 307, 307, 307, 307, 307, 307, 308, 308, 308, 308, 308, 308, 308, 309, 309, 310, 310, 310, 311, 311, 311, 311, 311,
            311, 312, 312, 312, 312, 312, 313, 313, 313, 313, 313, 314, 314, 314, 314, 314, 315, 315, 315, 315, 315, 315, 316, 316, 316, 316, 316,
            317, 317, 317, 317, 317, 317, 317, 318, 318, 318, 318, 319, 319, 319, 319, 319, 319, 319, 319, 319, 319, 320, 320, 320, 320, 320, 320,
            320, 320, 320, 321, 321, 321, 321, 321, 321, 321, 322, 322, 322, 322, 322, 322, 323, 323, 323, 323, 323, 324, 324, 324, 324, 325, 325,
            325, 325, 326, 327, 327, 328, 328, 329, 329, 329, 329, 330, 330, 330, 330, 330, 330, 330, 330, 331, 331, 331, 331, 331, 331, 332, 332,
            332, 332, 333, 333, 333, 333, 333, 333, 333, 333, 333, 334, 334, 334, 334, 335, 335, 335, 335, 335, 336, 336, 336, 336, 336, 337, 337,
            337, 337, 337, 337, 338, 338, 338, 338, 338, 338, 338, 339, 339, 340, 340, 340, 341, 341, 341, 341, 341, 342, 342, 342, 342, 342, 342,
            342, 342, 343, 344, 344, 344, 344, 344, 344, 345, 345, 346, 346, 346, 346, 346, 347, 347, 347, 347, 348, 348, 348, 348, 348, 348, 348,
            348, 349, 349, 349, 349, 350, 350, 350, 350, 350, 351, 351, 352, 352, 352, 352, 354, 354, 354, 354, 354, 355, 355, 355, 355, 356, 356,
            356, 356, 356, 357, 357, 358, 358, 359, 359, 359, 359, 361, 361, 361, 361, 362, 362, 362, 363, 363, 363, 363, 364, 364, 364, 364, 365,
            365, 365, 365, 365, 365, 366, 366, 366, 366, 367, 367, 367, 369, 369, 370, 370, 370, 370, 370, 371, 371, 372, 372, 372, 372, 372, 372,
            374, 374, 374, 374, 375, 376, 376, 376, 377, 377, 377, 378, 378, 379, 379, 379, 379, 380, 380, 380, 380, 381, 381, 382, 382, 383, 384,
            384, 385, 385, 385, 385, 386, 387, 387, 387, 388, 388, 389, 389, 389, 390, 390, 390, 391, 391, 392, 393, 394, 394, 396, 396, 396, 397,
            397, 397, 398, 398, 398, 399, 399, 399, 399, 399, 399, 400, 400, 402, 402, 402, 403, 403, 403, 403, 404, 404, 404, 404, 406, 407, 407,
            408, 408, 408, 409, 409, 409, 410, 410, 411, 411, 412, 412, 414, 416, 416, 417, 417, 417, 417, 419, 420, 420, 420, 420, 421, 421, 421,
            422, 422, 423, 423, 423, 423, 424, 424, 424, 425, 425, 425, 427, 427, 429, 429, 430, 430, 430, 430, 431, 433, 434, 435, 435, 435, 436,
            437, 437, 440, 440, 441, 442, 443, 444, 445, 446, 447, 448, 450, 450, 450, 450, 451, 452, 454, 455, 455, 455, 457, 458, 459, 459, 461,
            461, 463, 465, 466, 466, 467, 470, 471, 473, 474, 479, 479, 480, 481, 482, 483, 485, 491, 493, 493, 495, 496, 497, 498, 501, 507, 511,
            512, 512, 512, 515, 519, 521, 521, 521, 523, 524, 527, 531, 535, 537, 541, 542, 550, 556, 557, 557, 559, 561, 569, 572, 573, 575, 579,
            580, 580, 585, 589, 622, 627, 647, 667, 715, 757, 780, 1302, 4224, 10520, 16912, 17605, 17773, 18141, 18234, 19908, 22806, 26071, 
            33167, 39460, 39642, 55666, 59773, 63297, 74729, 82293, 84261, 103051, 125968, 126638, 133620, 134000, 156120, 156834, 158566, 159945,
            163617, 168131, 170576, 186151, 187063, 196981, 199264, 205006, 211618, 214498, 232698, 241394, 242280, 253799, 254397, 257196, 261598,
            316979, 329457, 654834};


        //Make an empty minimizer mapper just to get the fragment length distr
        gbwtgraph::GBWTGraph gbwt_graph;
        vector<gbwtgraph::DefaultMinimizerIndex*> minimizer_indexes;
        MinimumDistanceIndex distance_index;
        PathPositionHandleGraph* handle_graph;
        TestMinimizerMapper test_mapper (gbwt_graph, minimizer_indexes, distance_index, handle_graph);
        for (int64_t dist : distances) {
            test_mapper.fragment_length_distr.register_fragment_length(dist);
        }

        SECTION("Make sure we have a reasonable distribution") {
            //Distribution should ignore outliers
            REQUIRE(test_mapper.fragment_length_distr.std_dev() <= 400);
        }
}





}

}

