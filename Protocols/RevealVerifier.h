#ifndef PROTOCOLS_REVEAL_VERIFIER_H_
#define PROTOCOLS_REVEAL_VERIFIER_H_

#include "Math/gfpScalar.h"
#include "Math/gfpMatrix.h"

namespace hmmpc
{
// Batch check the correctness of the revealed values with t-sharing shares. 
class RevealVerifier
{
public:
    // Input
    static size_t reveal_num;
    static vector<gfpScalar> opened;
    static vector<gfpScalar> shared;
    static vector<gfpMatrix> openedMatrix;
    static vector<gfpMatrix> sharedMatrix;

    // Set
    static void push_scalar(gfpScalar openValue, gfpScalar shareValue);
    static void push_matrix(gfpMatrix openValues, gfpMatrix shareValues);

    // Verify
    static void batch_reveal_check();

    // Tool
    static gfpMatrix vector_2_matrix(vector<gfpScalar> v);
};

}

#endif