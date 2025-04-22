#pragma once
#include "Types/maliciousCheck.h"
#include "Tools/Subroutines.h"
#include "Protocols/ShareBundle.h"
#include "Protocols/RandomShare.h"
using Eigen::seq;
using Eigen::seqN;
using Eigen::last;

namespace hmmpc
{

// Falcon
void funcCheckMaliciousMatMulFalcon(const gfpMatrix &x, const gfpMatrix &y, const gfpMatrix &z)
{
    size_t rows = x.rows();
    size_t common_dim = x.cols();
    size_t columns = y.cols();
    ShareBundle tmp(x.size()+y.size(), 1);
    BeaverBundleTriple trip(rows, common_dim, columns);
    trip.beaver_triple();

    tmp.shares(seqN(0, x.size()), 0) = (x - trip.left_shares).reshaped<RowMajor>();
    tmp.shares(seqN(x.size(), y.size()), 0) = (y - trip.right_shares).reshaped<RowMajor>();

    tmp.reveal(true);

    gfpMatrix o(rows, columns);
    gfpMatrix u = tmp.secret()(seqN(0, x.size()), 0).reshaped<RowMajor>(rows, common_dim);
    gfpMatrix v = tmp.secret()(seqN(x.size(), y.size()), 0).reshaped<RowMajor>(common_dim, columns);
    o = trip.prod_shares + (u * trip.right_shares) + (trip.left_shares * v) + (u * v) - z;

    // Fcoin
    octet seed[SEED_SIZE];
    Create_Random_Seed(seed, (*ShareBase::P), SEED_SIZE);
    PRNG G;
    G.SetSeed(seed);
    gfpMatrix coefficients(rows*columns, 1);
    random_matrix(coefficients, G);
    
    Share opened;
    opened.share = (o.reshaped<RowMajor>().array() * coefficients.array()).sum();
    gfpScalar opendO = opened.reveal_with_check();
    if (opendO != gfpScalar(0)) {
        throw mult_verify_fail();
    }
}

void funcCheckMaliciousWiseMulFalcon(const gfpMatrix &x, const gfpMatrix &y, const gfpMatrix &z)
{
    size_t rows = x.rows();
    size_t columns = x.cols();
    ShareBundle tmp(x.size()<<1, 1);
    gfpMatrix A(rows, columns);
    gfpMatrix B(rows, columns);
    gfpMatrix C(rows, columns);
    A.setConstant(0);
    B.setConstant(0);
    C.setConstant(0); // C.array() = A.array() * B.array()

    tmp.shares(seqN(0, x.size()), 0) = (x.array() - A.array()).reshaped<RowMajor>();
    tmp.shares(seqN(x.size(), y.size()), 0) = (y.array() - B.array()).reshaped<RowMajor>();

    tmp.reveal(true);

    gfpMatrix o(rows, columns);
    gfpMatrix u = tmp.secret()(seqN(0, x.size()), 0).reshaped<RowMajor>(rows, columns);
    gfpMatrix v = tmp.secret()(seqN(x.size(), y.size()), 0).reshaped<RowMajor>(rows, columns);
    o = C.array() + (u.array() * B.array()) + (v.array() * A.array()) + (u.array() * v.array()) - z.array();

    // Fcoin
    octet seed[SEED_SIZE];
    Create_Random_Seed(seed, (*ShareBase::P), SEED_SIZE);
    PRNG G;
    G.SetSeed(seed);
    gfpMatrix coefficients(rows*columns, 1);
    random_matrix(coefficients, G);
    
    Share opened;
    opened.share = (o.reshaped<RowMajor>().array() * coefficients.array()).sum();
    gfpScalar opendO = opened.reveal_with_check();
    if (opendO != gfpScalar(0)) {
        throw mult_verify_fail();
    }
}

void funcCheckMaliciousMatMulLN17(const gfpMatrix &x, const gfpMatrix &y, const gfpMatrix &z)
{
    size_t rows = x.rows();
    size_t common_dim = x.cols();
    size_t columns = y.cols();
    ShareBundle tmp(x.size()+y.size(), 1);
    BeaverBundleTriple trip(rows, common_dim, columns);
    trip.beaver_triple();

    // Fcoin
    octet seed[SEED_SIZE];
    Create_Random_Seed(seed, (*ShareBase::P), SEED_SIZE);
    PRNG G;
    G.SetSeed(seed);
    gfpMatrix coefficients(rows*columns+1, 1);
    random_matrix(coefficients, G);
    gfpScalar alpha = coefficients(0, 0);

    tmp.shares(seqN(0, x.size()), 0) = (alpha * x.array() + trip.left_shares.array()).reshaped<RowMajor>();
    tmp.shares(seqN(x.size(), y.size()), 0) = (y + trip.right_shares).reshaped<RowMajor>();

    tmp.reveal(true);

    gfpMatrix o(rows, columns);
    gfpMatrix u = tmp.secret()(seqN(0, x.size()), 0).reshaped<RowMajor>(rows, common_dim);
    gfpMatrix v = tmp.secret()(seqN(x.size(), y.size()), 0).reshaped<RowMajor>(common_dim, columns);
    o = alpha * z.array() - trip.prod_shares.array() + (u * (trip.right_shares - v)).array() + (trip.left_shares * v).array();
    
    Share opened;
    opened.share = (o.reshaped<RowMajor>().array() * coefficients(seqN(1, rows*columns), 0).array()).sum();
    gfpScalar opendO = opened.reveal_with_check();
    if (opendO != gfpScalar(0)) {
        throw mult_verify_fail();
    }
}
}