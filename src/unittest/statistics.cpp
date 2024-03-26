/// \file statistics.cpp
///  
/// unit tests for statistics functions
///

#include "vg/io/json2pb.h"
#include "random_graph.hpp"
#include "catch.hpp"

#include "statistics.hpp"

#include <cmath>

namespace vg {
namespace unittest {

TEST_CASE("Summary statistics", "[statistics]") {
    SECTION("Even number of integers") {
        std::map<size_t, size_t> values { { 0, 1 }, { 1, 2 }, { 3, 1 }, { 4, 1 }, { 6, 1 } };
        SummaryStatistics statistics = summary_statistics(values);
        REQUIRE(statistics.mean == 2.5);
        REQUIRE(statistics.median == 2.0);
        REQUIRE(statistics.stdev == std::sqrt(25.5 / 6));
        REQUIRE(statistics.mode == 1.0);
        REQUIRE(statistics.number_of_values == 6);
        REQUIRE(statistics.max_value == 6.0);
        REQUIRE(statistics.count_of_max == 1);
    }

    SECTION("Odd number of doubles") {
        std::map<double, size_t> values { { 1.0, 2 }, { 3.0, 2 }, { 5.0, 3 }, { 10.0, 1 }, { 12.0, 1 } };
        SummaryStatistics statistics = summary_statistics(values);
        REQUIRE(statistics.mean == 5.0);
        REQUIRE(statistics.median == 5.0);
        REQUIRE(statistics.stdev == std::sqrt(114.0 / 9));
        REQUIRE(statistics.mode == 5.0);
        REQUIRE(statistics.number_of_values == 9);
        REQUIRE(statistics.max_value == 12.0);
        REQUIRE(statistics.count_of_max == 1);
    }
}

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

TEST_CASE("Max exponential fitting function learns parameters that are approximately correct", "[statistics]") {
    
    // data generated from a MaxExp(1/3, 5)
    vector<double> data{5.337616010050802, 3.915297109298315, 8.560560253325722, 3.442636743851616, 7.351206332454097, 3.2457226322649593, 4.764915780123888, 2.196821688140279, 14.671509268264463, 3.369084045294163, 7.539591129796632, 7.89165173016741, 3.758163876824844, 12.287663090879452, 8.695348720410442, 16.364009230530307, 2.7098247861759934, 5.254071473849081, 4.350829856734985, 8.469835362964952, 7.691679807303741, 3.823696015327002, 7.094947764999837, 2.7997745681160713, 7.482114390192067, 2.8239943130818337, 2.4486066178959636, 11.498672862221834, 12.337888440727223, 3.082956864545454, 8.37930971049281, 5.193272313161673, 5.484463085760223, 5.295411998740785, 6.868267520696233, 6.492116230197649, 3.362706049041917, 6.843369312683012, 5.487460294757351, 6.458508393744036, 4.907175735190073, 6.073701992306227, 9.716505070641006, 3.9371689355823447, 3.6440891662468227, 1.8084379689501324, 8.296064877871718, 6.675999080951752, 2.806886003823545, 4.060459885720452, 4.276134145787219, 12.143588685338079, 10.140203203267053, 17.129081471329716, 3.259610206716671, 6.764303746955554, 2.6119535388005017, 10.02253378177059, 11.841461956886059, 7.18105295271147, 10.901947016461397, 4.17064537010616, 5.522689522067273, 8.33790771032858, 10.391445792520498, 4.948468918372203, 11.498798449053126, 3.808757630459432, 6.326762840072416, 4.285478728667802, 5.177265464133489, 4.311588260899816, 2.306408237105522, 3.741723475740794, 6.428436394059439, 5.045538703872025, 3.543260332033181, 2.887563441308457, 3.6503129265076786, 6.749854406424635, 6.082987087420953, 4.861711767453076, 8.07095614547319, 9.760169342126193, 8.300926718401369, 5.5105564990107885, 3.0290928616036976, 1.9377822329492262, 2.1894612223848697, 6.464366187923735, 5.290498999138795, 3.3092341711092317, 9.021787003342192, 3.3963708298863824, 6.121722426380974, 9.144659299346355, 2.9098318763072584, 2.8694629720489386, 8.513678309706746, 8.803540814918778};
    
    auto params = fit_max_exponential(data);
    
    // more or less testing "does it look right" because of the very slack bounds
    
    REQUIRE(abs(params.first - .3333) < 1.0);
    REQUIRE(abs(params.second - 5.0) < 1.5);
}

TEST_CASE("Phred probability addition works", "[statistics]") {

    vector<double> phred1 = {prob_to_phred(0.01), prob_to_phred(1E-8), prob_to_phred(0.1), 30 };
    vector<double> phred2 = {prob_to_phred(0.1), prob_to_phred(0.3), prob_to_phred(0.2), 30 };
    vector<double> result = {prob_to_phred(0.11), prob_to_phred(0.3+1E-8), prob_to_phred(0.3), 26.98970004336};
    
    for (size_t i = 0; i < phred1.size(); i++) {
        REQUIRE((size_t) round(phred_add(phred1[i], phred2[i])) == (size_t) round(result[i]));
    }
    
}

TEST_CASE("Phred probability summation works", "[statistics]") {

    vector<double> to_sum = {30, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
    double result = 25.2287874528;
   
    REQUIRE(phred_sum(to_sum) == Approx(result));
    
}

TEST_CASE("Truncated normal functions produce correct results", "[statistics]") {
    
    random_device rd;
    default_random_engine prng(rd());
    
    for (int s = 1; s <= 3; ++s) {
        for (int m = -2; m <= 2; ++m) {
            for (int x = m-4; x <= m+4; ++x) {
                
                double d_norm = normal_pdf<double>(x, m, s);
                
                truncated_normal_distribution<> trunc1(m, s, -numeric_limits<double>::max(), m);
                truncated_normal_distribution<> trunc2(m, s, m, numeric_limits<double>::max());
                
                double d_trunc1 = trunc1.density(x);
                double d_trunc2 = trunc2.density(x);
                
                double p_trunc1 = trunc1.cumul(x);
                double p_trunc2 = trunc2.cumul(x);
                if (x < m) {
                    REQUIRE(abs((d_trunc1 / 2.0 - d_norm) / d_norm) < .0001);
                    REQUIRE(d_trunc2 == 0.0);
                    double p_norm1 = Phi(double(x - m) / s);
                    REQUIRE(abs(p_trunc1 / 2.0 - p_norm1) / p_norm1 < .0001);
                    REQUIRE(p_trunc2 == 0.0);
                }
                else if (x > m) {
                    REQUIRE(abs((d_trunc2 / 2.0 - d_norm) / d_norm) < .0001);
                    REQUIRE(d_trunc1 == 0.0);
                    double p_norm2 = Phi(double(x - m) / s) - .5;
                    REQUIRE(abs(p_trunc2 / 2.0 - p_norm2) / p_norm2 < .0001);
                    REQUIRE(p_trunc1 == 1.0);
                }
                else {
                    REQUIRE(abs(p_trunc1 - 1.0) < .0001);
                    REQUIRE(p_trunc2 < .0001);
                }
                
                for (int k = 0; k < 20; ++k) {
                    REQUIRE(trunc1(prng) <= m);
                    REQUIRE(trunc2(prng) >= m);
                }
            }
        }
    }
    
    for (int s = 1; s <= 3; ++s) {
        for (int m = -2; m <= 2; ++m) {
            for (int lo = m - 2; lo < m + 2; ++lo) {
                for (int hi = lo + 1; hi <= m + 2; ++hi) {
                    truncated_normal_distribution<> distr(m, s, lo, hi);
                    //cerr << "m " << m << " s " << s << " l " << lo << " h " << hi << " mean " << distr.mean() << endl;
                    REQUIRE(distr.stddev() < s);
                    if (hi <= m) {
                        REQUIRE(distr.mean() < m);
                    }
                    else if (lo >= m) {
                        REQUIRE(distr.mean() > m);
                    }
                    else if ((m - lo) == (hi - m)) {
                        REQUIRE(abs(distr.mean() - m) < .0001);
                    }
                    else if ((m - lo) > (hi - m)) {
                        REQUIRE(distr.mean() < m);
                    }
                    else {
                        REQUIRE(distr.mean() > m);
                    }
                    
                    for (int k = 0; k < 20; ++k) {
                        auto x = distr(prng);
                        REQUIRE(x <= hi);
                        REQUIRE(x >= lo);
                    }
                }
            }
        }
    }
}
}
}
