#include "Protocols/RevealVerifier.h"
#include "Protocols/RandomShare.h"

#include "Tools/Subroutines.h"
using Eigen::seq;
using Eigen::seqN;
using Eigen::last;

namespace hmmpc
{

// Initialize
size_t RevealVerifier::reveal_num = 0;
vector<gfpScalar> RevealVerifier::opened;
vector<gfpScalar> RevealVerifier::shared;
vector<gfpMatrix> RevealVerifier::openedMatrix;
vector<gfpMatrix> RevealVerifier::sharedMatrix;

void RevealVerifier::push_scalar(gfpScalar openValue, gfpScalar shareValue)
{
    opened.push_back(openValue);
    shared.push_back(shareValue);
    reveal_num++;
}

void RevealVerifier::push_matrix(gfpMatrix openValues, gfpMatrix shareValues)
{
    assert(openValues.cols() == shareValues.cols());
    assert(openValues.rows() == shareValues.rows());
    openedMatrix.push_back(openValues.reshaped<RowMajor>());
    sharedMatrix.push_back(shareValues.reshaped<RowMajor>());
    reveal_num += openValues.size();
}

// TODO_ZYS: May be adding the mutex lock when clear the opened and shared vectors
void RevealVerifier::batch_reveal_check()
{
    if (reveal_num <= 0) {
        return;
    }
    // Fcoin
    octet seed[SEED_SIZE];
    Create_Random_Seed(seed, (*ShareBase::P), SEED_SIZE);
    PRNG G;
    G.SetSeed(seed);
    gfpMatrix coefficients(reveal_num, 1);
    random_matrix(coefficients, G);

    // Combine
    gfpMatrix opened_matrix(reveal_num, 1);
    gfpMatrix shared_matrix(reveal_num, 1);
    opened_matrix(seqN(0, opened.size()), 0) = vector_2_matrix(opened);
    shared_matrix(seqN(0, opened.size()), 0) = vector_2_matrix(shared);
    size_t index = opened.size();
    for(size_t i = 0; i < openedMatrix.size(); i++) {
        opened_matrix(seqN(index, openedMatrix[i].size()), 0) = openedMatrix[i];
        shared_matrix(seqN(index, openedMatrix[i].size()), 0) = sharedMatrix[i];
        index += openedMatrix[i].size();
    }
    Share verified;
    verified.set_share((coefficients.array()*shared_matrix.array()).sum()-(coefficients.array()*opened_matrix.array()).sum());
    gfpScalar z = verified.reveal_with_check();

    // Check whether z is equal to zero
    if (z != gfpScalar(0)) {
        throw reveal_verify_fail();
    }

    // Clear
    opened.clear();
    shared.clear();
}

gfpMatrix RevealVerifier::vector_2_matrix(vector<gfpScalar> v)
{
    int num = v.size();
    gfpMatrix res(num, 1);
    for (int i = 0; i < num; i++)
    {
        res(i) = v[i];
    }
    return res;
}
}