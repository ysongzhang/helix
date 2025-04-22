#ifndef PROTOCOLS_INPUT_VERIFIER_H_
#define PROTOCOLS_INPUT_VERIFIER_H_

#include "Math/gfpScalar.h"
#include "Math/gfpMatrix.h"

namespace hmmpc
{

// Batch check the correctness of all the Input sharings in the online phase and all the t-sharings in the offline phase.
class InputVerifier
{
public:
    // Input
    static size_t share_num;
    static size_t share_num_offline;
    static vector<gfpScalar> shares;
    static vector<gfpScalar> shares_offline;
    static vector<gfpMatrix> shareMatrix;
    static vector<gfpMatrix> shareMatrix_offline;

    // Set
    static void push_scalar(gfpScalar share_scalar, bool is_offline = false);
    static void push_matrix(gfpMatrix share_matrix, bool is_offline = false);

    // Verify
    static void batch_correctness_check_on_shares(bool is_offline = false);

    // Tool
    static gfpMatrix vector_2_matrix(vector<gfpScalar> v);
};

}

#endif