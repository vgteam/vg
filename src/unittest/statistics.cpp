/// \file statistics.cpp
///  
/// unit tests for statistics functions
///

#include "json2pb.h"
#include "random_graph.hpp"
#include "catch.hpp"

#include "statistics.hpp"

#include <cmath>

namespace vg {
namespace unittest {

TEST_CASE("Matrix algebra functions work correctly", "[matrix][statistics]") {
    
    SECTION("Transpose produces correct output") {
        vector<vector<double>> mat{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
        auto mat_T = transpose(mat);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                REQUIRE(mat_T[i][j] == mat[j][i]);
            }
        }
    }
    
    SECTION("Matrix multiply produces correct output") {
        vector<vector<double>> A{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        vector<vector<double>> B{{1.0, 2.0}, {4.0, 5.0}, {7.0, 8.0}};
        vector<double> b{1.0, 2.0, 3.0};
        
        auto AB = matrix_multiply(A, B);
        auto Ab = matrix_multiply(A, b);
        REQUIRE(AB.size() == 3);
        REQUIRE(AB.front().size() == 2);
        REQUIRE(Ab.size() == 3);
        for (size_t i = 0; i < 3; ++i) {
            REQUIRE(Ab[i] == b[i]);
            for (size_t j = 0; j < 2; ++j) {
                REQUIRE(AB[i][j] == B[i][j]);
            }
        }
    }
    
    SECTION("Matrix inversion produces correct output on an identity matrix") {
        
        vector<vector<double>> A{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        
        auto A_inv = matrix_invert(A);
        auto B = matrix_multiply(A, A_inv);
        REQUIRE(B.size() == 3);
        REQUIRE(B.front().size() == 3);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                if (i == j) {
                    REQUIRE(abs(B[i][j] - 1.0) < 1e-6);
                }
                else {
                    REQUIRE(abs(B[i][j] - 0.0) < 1e-6);
                }
            }
        }
    }
    
    SECTION("Matrix inversion produces correct output on a dense matrix") {
                    
        vector<vector<double>> A{{1.0, 2.0, 3.0}, {4.0, 4.0, 6.0}, {7.0, 8.0, 9.0}};
        
        auto A_inv = matrix_invert(A);
        auto B = matrix_multiply(A, A_inv);
        REQUIRE(B.size() == 3);
        REQUIRE(B.front().size() == 3);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                if (i == j) {
                    REQUIRE(abs(B[i][j] - 1.0) < 1e-6);
                }
                else {
                    REQUIRE(abs(B[i][j] - 0.0) < 1e-6);
                }
            }
        }
    }
    
    SECTION("Matrix inversion produces correct output on matrix with 0 diagonal elements") {
        
        vector<vector<double>> A{{0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}};
        
        auto A_inv = matrix_invert(A);
        auto B = matrix_multiply(A, A_inv);
        REQUIRE(B.size() == 3);
        REQUIRE(B.front().size() == 3);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                if (i == j) {
                    REQUIRE(abs(B[i][j] - 1.0) < 1e-6);
                }
                else {
                    REQUIRE(abs(B[i][j] - 0.0) < 1e-6);
                }
            }
        }
    }
    
    SECTION("Regression produces correct output") {
        
        vector<vector<double>> X{{2.0, 4.0}, {1.0, 3.0}, {9.0, -1.0}};
        vector<double> y{-2.0, -2.0, 10.0};
        
        auto b = regress(X, y);
        REQUIRE(b.size() == 2);
        REQUIRE(abs(b[0] - 1.0) < 1e-6);
        REQUIRE(abs(b[1] + 1.0) < 1e-6);
    }
}

TEST_CASE("Golden section search finds the maximum value of functions", "[statistics]") {
    
    double tol = 1e-7;
    
    double x_star = golden_section_search([](double x){return -x * x;}, -5.0, 3.0, tol);
    REQUIRE(abs(x_star - 0.0) < tol);
    x_star = golden_section_search([](double x){return -(x - 1.0) * (x - 1.0);}, -5.0, 3.0, tol);
    REQUIRE(abs(x_star - 1.0) < tol);
}

TEST_CASE("Weibull fitting function learns parameters that are approximately correct", "[statistics]") {
    
    // data generated from a Weibull(2, 5, -5)
    vector<double> data{-3.2863010560586057, -3.265901197299282, -3.2395977179010704, -3.298833227141367, -2.615113915160391, -3.9324298194158724, -4.024530179324111, -2.59497751154026, -3.2548199030651226, -2.9467899257630905, -2.9799180391600166, -3.4749486393864917, -3.7315849540243677, -2.63469189535821, -2.8219021918011764, -3.358173575916232, -2.9383155088608266, -3.379009701282137, -3.260318275682521, -2.9490579397824965, -3.4351178840323624, -2.602171270366019, -2.803953248512817, -4.065155173527343, -3.2839353893572643, -2.543548895999789, -2.5121960265618855, -2.38467202060086, -3.0120113277973433, -3.100453455801046, -3.2085737188843613, -2.734073830714494, -3.9077865209223948, -3.6713719145120707, -3.422549578221937, -3.5246480818245356, -3.141009851431637, -3.2855874557475557, -2.773266164013716, -3.2543719223151495, -4.269537957765088, -3.4469570357185546, -2.9427257648842873, -3.7528826052108064, -3.3325105543515727, -3.234106579900243, -3.2565518957851305, -3.926859392748837, -2.8189608032656333, -3.028763054440476, -3.4099551565220567, -3.431991762459708, -3.4278906383226175, -3.08233531351979, -3.352921834952518, -3.1324533067172613, -3.021674739771419, -2.6581471398649232, -3.685320035803991, -3.836252464389406, -2.9593194188413574, -3.393645770688131, -3.3859553717377118, -2.8574196914091394, -2.321850957638378, -3.4047506312834592, -2.1433410782306663, -3.587814941158575, -3.607956803220886, -2.9343208317863674, -3.153120697794097, -3.2981588087237528, -2.9059552752804922, -3.7653847530862166, -2.4576008046142968, -2.7262542404702823, -2.7428319651236364, -3.207397194907559, -2.5086102154630985, -2.8838124024584135, -3.4714216443316372, -3.308274277914685, -2.624649189844103, -3.3206285066384686, -2.9895455335460155, -2.217439537389177, -2.9889485346452975, -3.82827474544077, -3.725429515982807, -2.671231052001676, -3.497747321877804, -3.580491178009275, -2.5262955303405272, -3.0248001671970526, -3.2195595191796764, -2.558484876025462, -2.599339238080604, -2.7237085963716945, -3.306923928993285, -3.822948500388329};
    
    auto params = fit_offset_weibull(data);
    
    // more or less testing "does it look right" because of the very slack bounds
    
    //cerr << get<0>(params) << " " << get<1>(params) << " " << get<2>(params) << endl;
    REQUIRE(abs(get<0>(params) - 2.0) < 2.0);
    REQUIRE(abs(get<1>(params) - 5.0) < 2.0);
    REQUIRE(abs(get<2>(params) + 5.0) < 2.0);
}


}
}
