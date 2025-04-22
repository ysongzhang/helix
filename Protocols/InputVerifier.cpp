#include "Protocols/InputVerifier.h"
#include "Protocols/RandomShare.h"

#include "Tools/Subroutines.h"
using Eigen::seq;
using Eigen::seqN;
using Eigen::last;

namespace hmmpc
{

// Initialize
size_t InputVerifier::share_num = 0;
size_t InputVerifier::share_num_offline = 0;
vector<gfpScalar> InputVerifier::shares;
vector<gfpScalar> InputVerifier::shares_offline;
vector<gfpMatrix> InputVerifier::shareMatrix;
vector<gfpMatrix> InputVerifier::shareMatrix_offline;

void InputVerifier::push_scalar(gfpScalar share_scalar, bool is_offline)
{
    if(is_offline) {
        shares_offline.push_back(share_scalar);
        share_num_offline++;
    } else{
        shares.push_back(share_scalar);
        share_num++;
    }
}

void InputVerifier::push_matrix(gfpMatrix share_matrix, bool is_offline)
{
    if(is_offline) {
        shareMatrix_offline.push_back(share_matrix.reshaped<RowMajor>());
        share_num_offline += share_matrix.size();
    } else{
        shareMatrix.push_back(share_matrix.reshaped<RowMajor>());
        share_num += share_matrix.size();
    }
}

// TODO_ZYS: May be adding the mutex lock when clear the shares vector
void InputVerifier::batch_correctness_check_on_shares(bool is_offline)
{
    if(is_offline){
        if (share_num_offline <= 0) {
            return;
        }

        // Fcoin
        octet seed[SEED_SIZE];
        Create_Random_Seed(seed, (*ShareBase::P), SEED_SIZE);
        PRNG G;
        G.SetSeed(seed);
        gfpMatrix coefficients(share_num_offline, 1);
        random_matrix(coefficients, G);

        // Frand
        gfpMatrix r(1, 1);
        RandomShare::generate_passive_random_sharings_PRG(r);

        // Combine
        gfpMatrix shared_matrix(share_num_offline, 1);
        shared_matrix(seqN(0, shares_offline.size()), 0) = vector_2_matrix(shares_offline);
        size_t index = shares_offline.size();
        for(size_t i = 0; i < shareMatrix_offline.size(); i++) {
            shared_matrix(seqN(index, shareMatrix_offline[i].size()), 0) = shareMatrix_offline[i];
            index += shareMatrix_offline[i].size();
        }
        Share verified;
        verified.set_share((coefficients.array()*shared_matrix.array()).sum() + r(0, 0));
        verified.reveal_with_check();

        // Clear
        shares_offline.clear();
    } else{
        if (share_num <= 0) {
            return;
        }

        // Fcoin
        octet seed[SEED_SIZE];
        Create_Random_Seed(seed, (*ShareBase::P), SEED_SIZE);
        PRNG G;
        G.SetSeed(seed);
        gfpMatrix coefficients(share_num, 1);
        random_matrix(coefficients, G);

        // Frand
        gfpMatrix r(1, 1);
        RandomShare::generate_passive_random_sharings_PRG(r);

        // Combine
        gfpMatrix shared_matrix(share_num, 1);
        shared_matrix(seqN(0, shares.size()), 0) = vector_2_matrix(shares);
        size_t index = shares.size();
        for(size_t i = 0; i < shareMatrix.size(); i++) {
            shared_matrix(seqN(index, shareMatrix[i].size()), 0) = shareMatrix[i];
            index += shareMatrix[i].size();
        }
        Share verified;
        verified.set_share((coefficients.array()*shared_matrix.array()).sum() + r(0, 0));
        verified.reveal_with_check();

        // Clear
        shares.clear();
    }
}

gfpMatrix InputVerifier::vector_2_matrix(vector<gfpScalar> v) {
    int num = v.size();
    gfpMatrix res(num, 1);
    for (int i = 0; i < num; i++)
    {
        res(i) = v[i];
    }
    return res;
}

}