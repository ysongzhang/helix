#ifndef PROTOCOLS_MULT_VERIFIER_H_
#define PROTOCOLS_MULT_VERIFIER_H_

#include "Types/Triple.h"

namespace hmmpc
{

class MultVerifier;

const int COMPRESSION_PARAMETER = 5;

class MultVerifier
{
public:    
    // Type
    static string verify_protocol;
    static bool offline_flag;

    /// Input
    static size_t maxRandomNum;  // Store the max number of reused factor.
    static size_t maxRandomNum_offline;
    // Each element (gfpMatrix) contains 3 columns and each row represents a beaver triple.
    static vector<gfpMatrix> vectorTriples; 
    // Each element (gfpMatrix) contains n columns that represent n dot triples. 
    // For a column vector, the first (rows()-1)/2 elements are left vector, the second (rows()-1)/2 elements are right vector, and the last element is the product.
    static vector<gfpMatrix> vectorDotTriples; 
    // Each element (gfpMatrix) contains 2*m+1 columns and m represents the reused factor.
    // Each row represents a set of reused triples.
    static vector<gfpMatrix> vectorReusedTriples;
    // Each element (MatrixTriple) represents a matrix triple.
    static vector<MatrixTriple> vectorMatrixTriples;
    // Used in the offline phase
    static vector<gfpMatrix> vectorTriples_offline; 
    // Used in the offline phase
    static vector<gfpMatrix> vectorDotTriples_offline;
    // Used in the offline phase
    static vector<gfpMatrix> vectorReusedTriples_offline;
    // Used in the offline phase
    static vector<MatrixTriple> vectorMatrixTriples_offline;

    // Recurive
    static int k;
    static int triples_num;
    static int triples_num_offline;
    static gfpVector recurive_left;
    static gfpVector recurive_right;
    static gfpScalar recurive_product;

    static gfpMatrix reconstruction_matrix; // Use in the step two of verify

    // Set
    static void push_triple(const gfpScalar a, const gfpScalar b, const gfpScalar c, bool is_offline = false);
    static void push_triples(const gfpMatrix a, const gfpMatrix b, const gfpMatrix c, bool is_offline = false); // Note: The size of a, b and c must be equal.
    static void push_reused_triples(const gfpMatrix a, const gfpMatrix b, const gfpMatrix c, bool is_offline = false);
    static void push_dot_triple(const gfpVector a, const gfpVector b, const gfpScalar c, bool is_offline = false);
    static void push_matrix_triple(const gfpMatrix a, const gfpMatrix b, const gfpMatrix c, bool is_offline = false);

    // Verify
    static void mult_verify(bool is_offline = false);

    // Verify with [GS20]
    static void mult_verify_GS20();
    static void mult_verify_step_one();
    static void mult_verify_step_two();
    static void mult_verify_step_three();
    static gfpMatrix dimension_reduction(gfpMatrix input, size_t m);  // Those functions that return class may change to return void and use reference.
    static gfpMatrix extend_compress(gfpMatrix input, int flag, gfpScalar &reducedSum);
    static gfpMatrix compress(gfpMatrix input);
    static gfpMatrix combine_opt(gfpMatrix a, gfpMatrix b);  // Directly add.

    // Verify with optimized [LN17]
    static void mult_verify_LN17();
    static void mult_verify_with_triples();

    // Verify with the method in Falcon
    static void mult_verify_Falcon();
    static void mult_verify_with_correct_triples();

    // Init in set-up phase; Note that before use MultVerifier, must invoke the init()
    static void init();  // This is to generate reconstruct matrix with size k-1 * k
    static void init(int compression_parameter);  // This is to generate reconstruct matrix with size k-1 * k and set the k value.
    static void init(string protocol);
    static void init(int compression_parameter, string protocol);

    // Tools
    static void coin(gfpMatrix &res);
    static void random(gfpMatrix &res);
    static gfpMatrix transform(vector<gfpScalar> input);
    static gfpScalar dot_product(gfpMatrix left, gfpMatrix right);
    static void dot_product_opt(const gfpMatrix &left, const gfpMatrix &right, gfpMatrix &res);
    static void dot_product_opt(gfpMatrix &input, size_t colNum);
    static void init_coefficients(gfpVector &res, gfpScalar base_value, size_t pow_index);

    // Test
    static bool triples_num_compare_flag;
    static int before_triples_num;
    static int before_triples_num_offline;
    static void print_triples_num();
    static void test_verifier();
};

}

#endif