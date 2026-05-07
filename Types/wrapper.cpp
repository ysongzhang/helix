#pragma once
#include "Types/wrapper.h"
#include "Types/maliciousCheck.h"
#include "Protocols/MultVerifier.h"

namespace hmmpc
{
void funcMatMul(const sfixMatrix&a, const sfixMatrix&b, sfixMatrix &res, bool a_transpose, bool b_transpose, size_t precision)
{
    gfpMatrix left;
    gfpMatrix right;
    if(!a_transpose && !b_transpose){
        res.share() = a.share() * b.share();
        left = a.share();
        right = b.share();
    } else if(!a_transpose && b_transpose){
        res.share() = a.share() * b.share().transpose();
        left = a.share();
        right = b.share().transpose();
    } else if(a_transpose && !b_transpose){
        res.share() = a.share().transpose() * b.share();
        left = a.share().transpose();
        right = b.share();
    } else {
        res.share() = a.share().transpose() * b.share().transpose();
        left = a.share().transpose();
        right = b.share().transpose();
    }
    if(precision==FIXED_PRECISION)
        res.reduce_truncate();
    else
        res.reduce_truncate(precision);

    if(MultVerifier::verify_protocol.compare("GS20") == 0) {
        // Version 1: defer check and use [GS20] method.
        MultVerifier::push_matrix_triple(left, right, res.checked());
    } else if(MultVerifier::verify_protocol.compare("LN17") == 0) {
        // Version 2: use the method in LN17 to check the Mat MUL
        funcCheckMaliciousMatMulLN17(left, right, res.checked());
    } else if(MultVerifier::verify_protocol.compare("Falcon") == 0) {
        // Version 3: use the method in Falcon to check the Mat MUL
        funcCheckMaliciousMatMulFalcon(left, right, res.checked());
    } else {
    }

    // Version 4: use Beaver triplets to compute MatMul
    // if(!a_transpose && !b_transpose){
    //     BeaverBundleTriple trip(a.rows(), a.cols(), b.cols());
    //     trip.beaver_triple();
    //     res.share() = trip.beaver_mult(a.share(), b.share());
    // } else if(!a_transpose && b_transpose){
    //     BeaverBundleTriple trip(a.rows(), a.cols(), b.rows());
    //     trip.beaver_triple();
    //     res.share() = trip.beaver_mult(a.share(), b.share().transpose());
    // } else if(a_transpose && !b_transpose){
    //     BeaverBundleTriple trip(a.cols(), a.rows(), b.cols());
    //     trip.beaver_triple();
    //     res.share() = trip.beaver_mult(a.share().transpose(), b.share());
    // } else {
    //     BeaverBundleTriple trip(a.cols(), a.rows(), b.rows());
    //     trip.beaver_triple();
    //     res.share() = trip.beaver_mult(a.share().transpose(), b.share().transpose());
    // }
    // if(precision==FIXED_PRECISION)
    //     res.truncate();
    // else
    //     res.truncate(precision);
    
}

void funcCwiseMul(const sfixMatrix&a, const sintMatrix &b, sfixMatrix &res)
{
    res.share() = a.share().array() * b.share().array();
    res.reduce_degree();
    MultVerifier::push_triples(a.share(), b.share(), res.share());
}

void funcDivision(const sfixMatrix&a, const sfixMatrix &b, sfixMatrix &res)
{
    res = divideRowwise(a, b);
}

void funcTrunc(sfixMatrix &res, size_t precision)
{
    res.truncate(precision);
}

void funcReLU(const sfixMatrix &input, sintMatrix &reluPrime, sfixMatrix &activations)
{
    // input.ReLU(reluPrime, activations);
    input.ReLU_opt(reluPrime, activations);
}

void funcOnlyReLU(const sfixMatrix &input, sfixMatrix &activations)
{
    // input.ReLU(activations);
    input.ReLU_opt(activations);
}

void funcConvMatMul(const sfixMatrix &a, const sfixMatrix &b, const sfixMatrix &biases, sfixMatrix &res,
                     size_t B, size_t oh, size_t ow, size_t Dout)
{
    // a: (B*ow*oh, f*f*Din)
    // b: (f*f*Din, Dout)

    // tmp = B*ow*oh, Dout
    // REPORT_ZYS: Could be replaced by tmp_res = a * b;
    sfixMatrix tmp_res(a.share() * b.share());
    tmp_res.reduce_truncate();

    if(MultVerifier::verify_protocol.compare("GS20") == 0) {
        // Version 1: defer check and use [GS20] method.
        MultVerifier::push_matrix_triple(a.share(), b.share(), tmp_res.checked());
    } else if(MultVerifier::verify_protocol.compare("LN17") == 0) {
        // Version 2: use the method in LN17 to check the Mat MUL
        funcCheckMaliciousMatMulLN17(a.share(), b.share(), tmp_res.checked());
    } else if(MultVerifier::verify_protocol.compare("Falcon") == 0) {
        // Version 3: use the method in Falcon to check the Mat MUL
        funcCheckMaliciousMatMulFalcon(a.share(), b.share(), tmp_res.checked());
    } else {
    }
    
    gfpMatrix tmp = tmp_res.share().array() + gfpMatrix(biases.share().colwise().replicate(a.rows())).array();  // REPORT_ZYS: BUG? Need to truncate first and then add the bias.
    // res = B, ow*oh*Dout
    size_t nRow=ow*oh;
    for(size_t i = 0; i < B; i++){
        size_t startRow=i*nRow;
        // BUG LOG: The channel is stored sequentially in one row.
        res.share().row(i) = tmp.middleRows(startRow, nRow).reshaped().transpose();
    }
    // res.share() = ( (a.share() * b.share()).array() + gfpMatrix(biases.share().colwise().replicate(a.rows())).array() ).reshaped(B, oh*ow*Dout);//colmajor
}

void funcMaxpool(const sfixMatrix &input, sintMatrix &maxPrime, sfixMatrix &activations, 
                size_t B, size_t Din, size_t oh, size_t ow, size_t f)
{
    // input: (B*ow*oh*Din, f*f)
    // activations: (B, ow*oh*Din)
    // maxPrime: (B, ow*oh*Din*f*f)

    sfixMatrix tmpActivations(B*ow*oh*Din, 1);
    sintMatrix tmpPrime(B*ow*oh*Din, f*f);
    tmpPrime.share().setConstant(1);
    // input.Maxpool(tmpPrime, tmpActivations);
    input.Maxpool_opt(tmpPrime, tmpActivations);
    activations.share() = tmpActivations.share().reshaped(B, ow*oh*Din);
    maxPrime.share() = tmpPrime.share().reshaped<RowMajor>(B, ow*oh*Din*f*f);
}

void funcOnlyMaxpool(const sfixMatrix &input, sfixMatrix &activations,
                 size_t B, size_t Din, size_t oh, size_t ow)
{
    // input: (B*ow*oh*Din, f*f)
    // activations: (B, ow*oh*Din)
    sfixMatrix tmpActivations(B*ow*oh*Din, 1);
    // input.Maxpool(tmpActivations);
    input.Maxpool_opt(tmpActivations);
    activations.share() = tmpActivations.share().reshaped<RowMajor>(B, ow*oh*Din);
}
}
