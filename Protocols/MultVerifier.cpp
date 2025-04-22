#include "Protocols/MultVerifier.h"
#include "Protocols/Share.h"
#include "Protocols/ShareBundle.h"
#include "Protocols/RandomShare.h"
#include "Tools/Subroutines.h"
using Eigen::seq;
using Eigen::seqN;
using Eigen::last;

namespace hmmpc
{

// Initialize
vector<gfpMatrix> MultVerifier::vectorTriples;
vector<gfpMatrix> MultVerifier::vectorDotTriples;
vector<gfpMatrix> MultVerifier::vectorReusedTriples;
vector<MatrixTriple> MultVerifier::vectorMatrixTriples;
vector<gfpMatrix> MultVerifier::vectorTriples_offline;
vector<gfpMatrix> MultVerifier::vectorDotTriples_offline;
vector<gfpMatrix> MultVerifier::vectorReusedTriples_offline;
vector<MatrixTriple> MultVerifier::vectorMatrixTriples_offline;
size_t MultVerifier::maxRandomNum = 1;
size_t MultVerifier::maxRandomNum_offline = 1;
int MultVerifier::triples_num = 0;
int MultVerifier::triples_num_offline = 0;
int MultVerifier::before_triples_num = 0;
int MultVerifier::before_triples_num_offline = 0;
gfpVector MultVerifier::recurive_left;
gfpVector MultVerifier::recurive_right;
gfpScalar MultVerifier::recurive_product(0);
int MultVerifier::k; // NOTE_ZYS: k must be larger than 1.
gfpMatrix MultVerifier::reconstruction_matrix;
string MultVerifier::verify_protocol = "";
bool MultVerifier::offline_flag = false;
bool MultVerifier::triples_num_compare_flag = false;


void MultVerifier::init()
{
    init(COMPRESSION_PARAMETER);
    verify_protocol = "GS20";
}

void MultVerifier::init(int compression_parameter)
{   
    assert(compression_parameter>1);
    k = compression_parameter;
    reconstruction_matrix.resize(k-1, k);
    ShareBase::init_reuse_factor_extend_compress(k-1);

    // Init reconstruction_matrix: invoke get_reconstruction_vector_general with reuse_factor k-1 times
    for (int i = 0; i < k-1; i++)
    {
        gfpScalar point(i + k);
        reconstruction_matrix.row(i) = ShareBase::get_reconstruction_vector_general(point, k-1, ShareBase::reuse_factor_extend_compress).transpose();
    }
    verify_protocol = "GS20";
}

void MultVerifier::init(string protocol)
{
    init(COMPRESSION_PARAMETER);
    if(protocol.compare("GS20") == 0 || protocol.compare("LN17") == 0 || protocol.compare("Falcon") == 0) {
        verify_protocol = protocol;
    } else {
        throw Verify_Protocol_Error(protocol + "is not implemented!");
    }
}

void MultVerifier::init(int compression_parameter, string protocol)
{
    init(compression_parameter);
    if(protocol.compare("GS20") == 0 || protocol.compare("LN17") == 0 || protocol.compare("Falcon") == 0) {
        verify_protocol = protocol;
    } else {
        throw Verify_Protocol_Error(protocol + "is not implemented!");
    }
}

void MultVerifier::push_triple(const gfpScalar a, const gfpScalar b, const gfpScalar c, bool is_offline)
{
    gfpMatrix t(1, 3);
    t(0) = a;
    t(1) = b;
    t(2) = c;
    if(is_offline){
        vectorTriples_offline.push_back(t);
        triples_num_offline++;
        if(triples_num_compare_flag) {
            before_triples_num_offline++;
        }
    } else{
        vectorTriples.push_back(t);
        triples_num++;
        if(triples_num_compare_flag) {
            before_triples_num++;
        }
    }
}

void MultVerifier::push_triples(const gfpMatrix a, const gfpMatrix b, const gfpMatrix c, bool is_offline)
{
    assert(a.rows()==b.rows() && a.rows()==c.rows());
    assert(a.cols()==b.cols() && a.cols()==c.cols());
    gfpMatrix t(a.size(), 3);
    t.col(0) = a.reshaped<RowMajor>();
    t.col(1) = b.reshaped<RowMajor>();
    t.col(2) = c.reshaped<RowMajor>();
    if(is_offline){
        vectorTriples_offline.push_back(t);
        triples_num_offline += a.size();
        if(triples_num_compare_flag) {
            before_triples_num_offline += a.size();
        }
    } else{
        vectorTriples.push_back(t);
        triples_num += a.size();
        if(triples_num_compare_flag) {
            before_triples_num += a.size();
        }
    }
}

void MultVerifier::push_reused_triples(const gfpMatrix a, const gfpMatrix b, const gfpMatrix c, bool is_offline)
{
    size_t reusedFactor = b.cols();
    size_t len = a.rows();
    assert(a.cols() == 1);
    assert(c.cols() == reusedFactor);
    assert(b.rows() == len && c.rows() == len);
    gfpMatrix t(len, (reusedFactor<<1)+1);
    t << a, b, c;
    if(is_offline){
        vectorReusedTriples_offline.push_back(t);
        triples_num_offline += len;
        maxRandomNum_offline = (reusedFactor > maxRandomNum_offline) ? reusedFactor : maxRandomNum_offline;
        if(triples_num_compare_flag) {
            before_triples_num_offline += b.size();
        }
    } else{
        vectorReusedTriples.push_back(t);
        triples_num += len;
        maxRandomNum = (reusedFactor > maxRandomNum) ? reusedFactor : maxRandomNum;
        if(triples_num_compare_flag) {
            before_triples_num += b.size();
        }
    }
}

void MultVerifier::push_dot_triple(const gfpVector a, const gfpVector b, const gfpScalar c, bool is_offline)
{
    size_t len = a.size();
    gfpMatrix t((len<<1)+1, 1);
    t(seqN(0, len), 0) = a;
    t(seqN(len, len), 0) = b;
    t(len<<1, 0) = c;
    if(is_offline) {
        vectorDotTriples_offline.push_back(t);
        triples_num_offline += len;
        if(triples_num_compare_flag) {
            before_triples_num_offline += len;
        }
    } else{
        vectorDotTriples.push_back(t);
        triples_num += len;
        if(triples_num_compare_flag) {
            before_triples_num += len;
        }
    }
}

void MultVerifier::push_matrix_triple(const gfpMatrix a, const gfpMatrix b, const gfpMatrix c, bool is_offline)
{
    MatrixTriple matTriple(a, b, c);
    size_t reusedFactor = c.rows() + c.cols();

    size_t len = a.cols();  // The size of compressed vector.
    if(is_offline) {
        vectorMatrixTriples_offline.push_back(matTriple);
        triples_num_offline += len;
        maxRandomNum_offline = (reusedFactor > maxRandomNum_offline) ? reusedFactor : maxRandomNum_offline;
        if(triples_num_compare_flag) {
            before_triples_num_offline += c.size() * len;
        }
    } else {
        vectorMatrixTriples.push_back(matTriple);
        triples_num += len;
        maxRandomNum = (reusedFactor > maxRandomNum) ? reusedFactor : maxRandomNum;
        if(triples_num_compare_flag) {
            before_triples_num += c.size() * len;
        }
    }
}

void MultVerifier::mult_verify(bool is_offline)
{
    offline_flag = is_offline;
    if(verify_protocol.compare("GS20") == 0) {
        mult_verify_GS20();
    } else if(verify_protocol.compare("LN17") == 0) {
        mult_verify_LN17();
    } else if(verify_protocol.compare("Falcon") == 0) {
        mult_verify_Falcon();
    } else {
        throw Verify_Protocol_Error(verify_protocol + "is not implemented!");
    }
}

// TODO_ZYS: May be adding the mutex lock when clear the vectors and recurive matrix
void MultVerifier::mult_verify_GS20()
{
    if(offline_flag) {
        if(triples_num_offline < 1) {
            return;
        } else{
            if(triples_num_compare_flag) {
                cout<<"Before triples number (offline): "<<before_triples_num_offline<<endl;
                cout<<"Compressed triples number (offline): "<<triples_num_offline<<endl;
            }

            mult_verify_step_one();
            mult_verify_step_two();
            mult_verify_step_three();   

            // Clear
            vectorTriples_offline.clear();
            vectorDotTriples_offline.clear();
            vectorReusedTriples_offline.clear();
            vectorMatrixTriples_offline.clear();
            maxRandomNum_offline = 1;
            triples_num_offline = 0;
            recurive_left.resize(0);
            recurive_right.resize(0);
            recurive_product = 0;
        }
    } else{
        if(triples_num < 1) {
            return;
        } else{
            if(triples_num_compare_flag) {
                cout<<"Before triples number (online): "<<before_triples_num<<endl;
                cout<<"Compressed triples number (online): "<<triples_num<<endl;
            }

            mult_verify_step_one();
            mult_verify_step_two();
            mult_verify_step_three();

            // Clear
            vectorTriples.clear();
            vectorDotTriples.clear();
            vectorReusedTriples.clear();
            vectorMatrixTriples.clear();
            maxRandomNum = 1;
            triples_num = 0;
            recurive_left.resize(0);
            recurive_right.resize(0);
            recurive_product = 0;
        }
    }
}

void MultVerifier::mult_verify_step_one()
{
    // step 1: generate random field elements by PRF+seed
    octet seed[SEED_SIZE];
    Create_Random_Seed(seed, (*ShareBase::P), SEED_SIZE);
    PRNG G;
    G.SetSeed(seed);

    if(offline_flag) {
        recurive_left.resize(triples_num_offline, 1);
        recurive_right.resize(triples_num_offline, 1);

        // step 2: Matrix compression and Reused triples compression
        gfpVector randomness(maxRandomNum_offline);
        random_matrix(randomness, G);
        for (MatrixTriple t : vectorMatrixTriples_offline)
        {
            size_t len = t.second_dimension();
            gfpVector compressed_left(len);
            gfpVector compressed_right(len);
            gfpScalar compressed_prod;
            t.compress(randomness(seqN(0, t.first_dimension())), randomness(seqN(t.first_dimension(), t.third_dimension())), compressed_left, compressed_right, compressed_prod);
            gfpMatrix padding((len<<1)+1, 1);
            padding(seqN(0, len), 0) = compressed_left;
            padding(seqN(len, len), 0) = compressed_right;
            padding(len<<1, 0) = compressed_prod;
            vectorDotTriples_offline.push_back(padding);
        }
        for (gfpMatrix t : vectorReusedTriples_offline)
        {
            size_t reusedFactor = (t.cols() - 1) / 2;
            size_t len = t.rows();
            gfpMatrix padding(len, 3);
            padding.setZero();
            padding.col(0) = t.col(0);
            for(size_t i = 0; i < reusedFactor; i++) {
                padding.col(1) += randomness(i) * t.col(i+1);
                padding.col(2) += randomness(i) * t.col(i+1+reusedFactor);
            }
            vectorTriples_offline.push_back(padding);
        }
        
        // step 3: Normalization
        size_t index = 0;
        for(gfpMatrix t : vectorTriples_offline) {
            size_t len = t.rows();
            gfpVector coefficients(len);
            random_matrix(coefficients, G);
            recurive_left(seqN(index, len), 0) = (t.col(0).array() * coefficients.array());
            recurive_right(seqN(index, len), 0) = t.col(1);
            recurive_product += (t.col(2).array() * coefficients.array()).sum();
            index += len;
        }
        for(gfpMatrix t : vectorDotTriples_offline) {
            size_t num = t.cols();
            gfpVector coefficients(num);
            random_matrix(coefficients, G);
            for (size_t i = 0; i < num; i++)
            {
                size_t len = (t.rows() - 1) / 2;
                recurive_left(seqN(index, len), 0) =  (t.col(i)(seqN(0, len), 0).array() * coefficients(i));
                recurive_right(seqN(index, len), 0) = t.col(i)(seqN(len, len), 0);
                recurive_product += (t.col(i)(len<<1, 0) * coefficients(i));
                index += len;
            }
        }
    } else{

        recurive_left.resize(triples_num, 1);
        recurive_right.resize(triples_num, 1);

        // step 2: Matrix compression and Reused triples compression
        gfpVector randomness(maxRandomNum);
        random_matrix(randomness, G);
        for (MatrixTriple t : vectorMatrixTriples)
        {
            size_t len = t.second_dimension();
            gfpVector compressed_left(len);
            gfpVector compressed_right(len);
            gfpScalar compressed_prod;
            t.compress(randomness(seqN(0, t.first_dimension())), randomness(seqN(t.first_dimension(), t.third_dimension())), compressed_left, compressed_right, compressed_prod);
            gfpMatrix padding((len<<1)+1, 1);
            padding(seqN(0, len), 0) = compressed_left;
            padding(seqN(len, len), 0) = compressed_right;
            padding(len<<1, 0) = compressed_prod;
            vectorDotTriples.push_back(padding);
        }
        for (gfpMatrix t : vectorReusedTriples)
        {
            size_t reusedFactor = (t.cols() - 1) / 2;
            size_t len = t.rows();
            gfpMatrix padding(len, 3);
            padding.setZero();
            padding.col(0) = t.col(0);
            for(size_t i = 0; i < reusedFactor; i++) {
                padding.col(1) += randomness(i) * t.col(i+1);
                padding.col(2) += randomness(i) * t.col(i+1+reusedFactor);
            }
            vectorTriples.push_back(padding);
        }

        // step 3: Normalization
        int index = 0;
        for(gfpMatrix t : vectorTriples) {
            size_t len = t.rows();
            gfpVector coefficients(len);
            random_matrix(coefficients, G);
            recurive_left(seqN(index, len), 0) = (t.col(0).array() * coefficients.array());
            recurive_right(seqN(index, len), 0) = t.col(1);
            recurive_product += (t.col(2).array() * coefficients.array()).sum();
            index += len;
        }
        for(gfpMatrix t : vectorDotTriples) {
            size_t num = t.cols();
            gfpVector coefficients(num);
            random_matrix(coefficients, G);
            for (size_t i = 0; i < num; i++)
            {
                size_t len = (t.rows() - 1) / 2;
                recurive_left(seqN(index, len), 0) =  (t.col(i)(seqN(0, len), 0).array() * coefficients(i));
                recurive_right(seqN(index, len), 0) = t.col(i)(seqN(len, len), 0);
                recurive_product += (t.col(i)(len<<1, 0) * coefficients(i));
                index += len;
            }
        }
    }
}

void MultVerifier::mult_verify_step_two()
{
    size_t in_dimension = recurive_left.size();
    gfpMatrix input((in_dimension<<1)+1, 1);
    input(seqN(0, in_dimension), 0) = recurive_left;
    input(seqN(in_dimension, in_dimension), 0) = recurive_right;
    input(in_dimension<<1) = recurive_product; 
    gfpMatrix output = dimension_reduction(input, in_dimension);
    int out_dimension = (output.size() - 1) / 2;
    recurive_left = output(seqN(0, out_dimension), 0);
    recurive_right = output(seqN(out_dimension, out_dimension), 0);
    recurive_product = output(out_dimension<<1);
}

void MultVerifier::mult_verify_step_three()
{
    int m = recurive_left.size();
    gfpMatrix input(3, m+1);
    ShareBundle cShareBundle(m, 1); // The 2t-sharings need to be reduce degree

    // F_rand
    gfpMatrix a0b0(2, 1);
    random(a0b0);
    // Set 1st column triple.
    input(0, 0) = a0b0(0);
    input(1, 0) = a0b0(1);

    // Set 2nd ~ (m+1)th colums triples.
    input(0, seqN(1, m)) = recurive_left.transpose();
    input(1, seqN(1, m)) = recurive_right.transpose();

    // Compute the product and Set the 3rd row product; The 2t-sharings need to be reduce degree
    input(2, 0) = a0b0(0) * a0b0(1);
    input(2, seqN(1, m-1)) = (recurive_left(seqN(0, m-1), 0).array() * recurive_right(seqN(0, m-1), 0).array()).transpose();
    input(2, m) = recurive_product;  // need to reduce the input(2, seqN(1, m-1))
    
    gfpMatrix output = compress(input);

    // reconstruct output with check
    ShareBundle outputShare = ShareBundle(3);
    outputShare.shares(0) = output(0);
    outputShare.shares(1) = output(1);
    outputShare.shares(2) = output(2);
    outputShare.reveal_with_check();
    if (outputShare.secret()(0, 0) * outputShare.secret()(1, 0) != outputShare.secret()(2, 0)) {
        throw mult_verify_fail();
    }
}

gfpMatrix MultVerifier::dimension_reduction(gfpMatrix input, size_t m)
{
    // step 1: Compute the length of block and the remainder length.
    int l = m / k;
    if (l < 2) {
        return input;  // Only here update the static seed
    }
    int r = m - (l * k);

    int flag = 0;  // Define the signal of combine
    if (r != 0) {
        flag = 1;
    }

    // step 2: Build blocks: triples_blk + triples_remain
    // triples_blk each column is a dot product triple and triples_remain is a column vector.
    gfpMatrix triples_blk(2*l+1, k);  // length = k, each element (each column) contains 2*l+1 elements.
    gfpMatrix triples_remain;  // may be, using flag, len = r < k
    for (int i = 0; i < k; i++)  // i*l ~ (i+1)*(l-1)
    {
        triples_blk(seqN(0, l), i) = input(seqN(i*l, l), 0);
        triples_blk(seqN(l, l), i) = input(seqN(i*l+m, l), 0);
    }
    dot_product_opt(triples_blk, k+flag-1);

    if (flag) {
        triples_remain.resize(2*r+1, 1);
        triples_remain(seqN(0, r), 0) = input(seqN(k*l, r), 0);
        triples_remain(seqN(r, r), 0) = input(seqN(k*l+m, r), 0);
    } else {
        triples_blk(2*l, k-1) = input(input.size()-1);  // need to reduce triples_blk(2*l, seqN(0, k+flag-1).sum)
    }
    
    
    // step 3: Extend compress
    gfpScalar reducedSum(0);
    gfpMatrix res = extend_compress(triples_blk, flag, reducedSum);

    if (flag) {
        triples_remain(2*r) = input(input.size()-1) - reducedSum;
        res = combine_opt(res, triples_remain);
    }
    res = dimension_reduction(res, l + r);  // return must satisfy (len / k) < 2

    return res;
}

// input: each column is a dot product triple. Contain k dot product triples (k columns).
// output: return compressed column vector.
gfpMatrix MultVerifier::extend_compress(gfpMatrix input, int flag, gfpScalar &reducedSum)
{
    assert(input.cols() == k);

    int l = (input.rows() - 1) / 2;
    gfpMatrix fAndGPoints(k, 2*l);
    gfpMatrix hPoints(2*k-1, 1);

    fAndGPoints = input(seqN(0, 2*l), seq(0, last)).transpose();

    ShareBundle zShareBundle(2*k-2+flag, 1); // The 2t-sharings need to be reduce degree; (k+flag-1)+(k-1)
    zShareBundle.shares(seqN(0, k+flag-1), 0) = input(2*l, seqN(0, k+flag-1)).transpose();
    
    gfpMatrix fAndGShares = reconstruction_matrix * fAndGPoints;  // The max run-time cost in our MultVerifier! The overall number of multiplication is 2*N*(k-1)
    zShareBundle.shares(seqN(k+flag-1, k-1), 0) = (fAndGShares.transpose()(seqN(0, l), seqN(0, k-1)).array() * 
                                                    fAndGShares.transpose()(seqN(l, l), seqN(0, k-1)).array())
                                                    .colwise().sum().transpose();

    zShareBundle.reduce_degree();  // reduce degree, only use semi-honest protocol

    // Fill hPoints
    reducedSum = zShareBundle.shares(seqN(0, k+flag-1), 0).sum();
    if(flag) {
        hPoints(seqN(0, 2*k-1), 0) = zShareBundle.shares;
    } else {
        hPoints(seqN(0, k-1), 0) = zShareBundle.shares(seqN(0, k-1), 0);
        hPoints(k-1, 0) = input(2*l, k-1) - reducedSum;
        hPoints(seqN(k, k-1), 0) = zShareBundle.shares(seqN(k-1, k-1), 0);
    }

    // Fcoin
    gfpMatrix tmp_r(2, 1);
    coin(tmp_r);
    gfpScalar r = tmp_r(0, 0);
    // assert(r > k-1); // Note_ZYS: The output of extend_compress is not opened. Thus, there is no need to abort when r <= k-1.

    // point == r
    gfpVector reconstructFAndG = ShareBase::get_reconstruction_vector_general(r, k-1, ShareBase::reuse_factor_extend_compress);
    gfpVector reconstructH = ShareBase::get_reconstruction_vector_general(r, 2*k-2);
    gfpMatrix fRAndGR = (reconstructFAndG.transpose() * fAndGPoints);  // size: 1 * 2l
    gfpMatrix hR = (reconstructH.transpose() * hPoints);  // size: 1 * 1

    gfpMatrix res(2*l+1, 1);  // column vector
    res(seqN(0, 2*l), 0) = fRAndGR.transpose();
    res(seqN(0, l), 0) = tmp_r(1, 0) * res(seqN(0, l), 0);
    res(2*l) = tmp_r(1, 0) * hR(0);
    return res;
}

// input: each column is a triple. Contain 3 rows.
// output: return compressed column vector.
gfpMatrix MultVerifier::compress(gfpMatrix input)
{
    int m = input.cols();  // the number of points
    gfpMatrix fAndGPoints(m, 2);
    gfpMatrix hPoints(2*m-1, 1);

    fAndGPoints = input(seqN(0, 2), seq(0, last)).transpose();
    
    ShareBundle zShareBundle(2*m-2, 1); // The 2t-sharings need to be reduce degree
    zShareBundle.shares(seqN(0, m-1), seq(0, last)) = input(2, seqN(0, m-1)).transpose();

    // reconstruct f(i)=x^i and h(i)=y^i, which are both degree-(m-1) polynomials.
    ShareBase::init_reuse_factor_compress(m-1);
    for(int point = m; point < (2 * m - 1); point++) {
        gfpVector reconstruct = ShareBase::get_reconstruction_vector_general(point, m-1, ShareBase::reuse_factor_compress);
        // reconstruct f(point) and h(point)
        gfpMatrix fAndGShare = (reconstruct.transpose() * fAndGPoints);  // size: 1 * 2
        zShareBundle.shares(point-1) = fAndGShare(0) * fAndGShare(1);
    }

    zShareBundle.reduce_degree();  // reduce degree, only use semi-honest protocol

    // Fill hPoints
    hPoints(seqN(0, m-1), 0) = zShareBundle.shares(seqN(0, m-1), seq(0, last));
    hPoints(m-1, 0) = input(2, m-1) - zShareBundle.shares(seqN(1, m-2), seq(0, last)).sum();
    hPoints(seqN(m, m-1), 0) = zShareBundle.shares(seqN(m-1, m-1), seq(0, last));

    // Fcoin
    gfpMatrix tmp_r(1, 1);
    coin(tmp_r);
    gfpScalar r = tmp_r(0, 0);
    assert(r > m-1);

    // point == r
    gfpVector reconstructFAndG = ShareBase::get_reconstruction_vector_general(r, m-1, ShareBase::reuse_factor_compress);
    gfpVector reconstructH = ShareBase::get_reconstruction_vector_general(r, 2*m-2);
    gfpMatrix fRAndGR = (reconstructFAndG.transpose() * fAndGPoints);
    gfpMatrix hR = (reconstructH.transpose() * hPoints);

    gfpMatrix res(3, 1);
    res(0) = fRAndGR(0);
    res(1) = fRAndGR(1);
    res(2) = hR(0);

    return res;
}

// Input: two column vectors, which means two dot product triples
// Output: one column vector
gfpMatrix MultVerifier::combine_opt(gfpMatrix a, gfpMatrix b)
{
    gfpMatrix res(a.size() + b.size() - 1, 1);
    int l_a = (a.size() - 1) / 2;
    int l_b = (b.size() - 1) / 2;
    
    res(seqN(0, l_a), 0) = a(seqN(0, l_a), 0);
    res(seqN(l_a, l_b), 0) = b(seqN(0, l_b), 0);
    res(seqN(l_a+l_b, l_a), 0) = a(seqN(l_a, l_a), 0);
    res(seqN(2*l_a+l_b, l_b), 0) = b(seqN(l_b, l_b), 0);
    res(2*(l_a+l_b)) = a(2*l_a) + b(2*l_b);
    return res;
}

void MultVerifier::mult_verify_LN17()
{
    if(offline_flag) {
        if(triples_num_offline < 1) {
            return;
        } else{
            mult_verify_step_one();
            mult_verify_with_triples();
            // Clear
            vectorTriples_offline.clear();
            vectorDotTriples_offline.clear();
            vectorReusedTriples_offline.clear();
            vectorMatrixTriples_offline.clear();
            maxRandomNum_offline = 1;
            triples_num_offline = 0;
            recurive_left.resize(0);
            recurive_right.resize(0);
            recurive_product = 0;
        }
    } else {
        if(triples_num < 1) {
            return;
        } else{

            mult_verify_step_one();
            mult_verify_with_triples();

            // Clear
            vectorTriples.clear();
            vectorDotTriples.clear();
            vectorReusedTriples.clear();
            vectorMatrixTriples.clear();
            maxRandomNum = 1;
            triples_num = 0;
            recurive_left.resize(0);
            recurive_right.resize(0);
            recurive_product = 0;
        }
        
    }
}

// For Test
void MultVerifier::mult_verify_with_triples()
{
    size_t in_dimension = recurive_left.size();
    ShareBundle tmp(2*in_dimension, 1);
    BeaverBundleTriple trip(1, in_dimension, 1);
    trip.beaver_triple(offline_flag);
    
    // Step 1: Fcoin
    gfpMatrix tmp_r(1, 1);
    coin(tmp_r);
    gfpScalar alpha = tmp_r(0, 0);

    // Step 2: Combine
    tmp.shares(seqN(0, in_dimension), 0) = alpha * recurive_left.array() + trip.left_shares.transpose().array();
    tmp.shares(seqN(in_dimension, in_dimension), 0) = recurive_right.array() + trip.right_shares.array();

    // Step 3: Open
    tmp.reveal_with_check();

    // Step 4: Compute and Open v
    Share v;
    v.share = alpha * recurive_product - trip.prod_shares(0, 0) + (tmp.secret()(seqN(0, in_dimension), 0).transpose() * trip.right_shares)(0, 0) + (trip.left_shares * tmp.secret()(seqN(in_dimension, in_dimension), 0))(0, 0) - (tmp.secret()(seqN(0, in_dimension), 0).transpose() * tmp.secret()(seqN(in_dimension, in_dimension), 0))(0, 0);
    gfpScalar opendV = v.reveal_with_check();
    if (opendV != gfpScalar(0)) {
        throw mult_verify_fail();
    }
}

void MultVerifier::mult_verify_Falcon()
{
    if(offline_flag) {
        if(triples_num_offline < 1) {
            return;
        } else{
            // Note_ZYS: Since the preprocessing phase of Falcon must provide the correct Beaver triples, we must use the GS20/LN17 protocol to verify the offline process.
            // Note_ZYS: Now use the GS20 protocol
            mult_verify_step_one();
            mult_verify_step_two();
            mult_verify_step_three();

            // Clear
            vectorTriples_offline.clear();
            vectorDotTriples_offline.clear();
            vectorReusedTriples_offline.clear();
            vectorMatrixTriples_offline.clear();
            maxRandomNum_offline = 1;
            triples_num_offline = 0;
            recurive_left.resize(0);
            recurive_right.resize(0);
            recurive_product = 0;
        }
    } else {
        if(triples_num < 1) {
            return;
        } else{
            mult_verify_step_one();
            mult_verify_with_correct_triples();

            // Clear
            vectorTriples.clear();
            vectorDotTriples.clear();
            vectorReusedTriples.clear();
            vectorMatrixTriples.clear();
            maxRandomNum = 1;
            triples_num = 0;
            recurive_left.resize(0);
            recurive_right.resize(0);
            recurive_product = 0;
        }
    }
}

void MultVerifier::mult_verify_with_correct_triples()
{
    size_t in_dimension = recurive_left.size();
    ShareBundle tmp(2*in_dimension, 1);
    BeaverBundleTriple trip(1, in_dimension, 1);
    trip.beaver_triple();

    tmp.shares(seqN(0, in_dimension), 0) = recurive_left.array() - trip.left_shares.transpose().array();
    tmp.shares(seqN(in_dimension, in_dimension), 0) = recurive_right.array() - trip.right_shares.array();
    tmp.reveal_with_check();
    Share v;
    v.share = trip.prod_shares(0, 0) + (tmp.secret()(seqN(0, in_dimension), 0).transpose() * trip.right_shares)(0, 0) + (trip.left_shares * tmp.secret()(seqN(in_dimension, in_dimension), 0))(0, 0) + (tmp.secret()(seqN(0, in_dimension), 0).transpose() * tmp.secret()(seqN(in_dimension, in_dimension), 0))(0, 0) - recurive_product;
    gfpScalar opendV = v.reveal_with_check();
    if (opendV != gfpScalar(0)) {
        throw mult_verify_fail();
    }
}

void MultVerifier::coin(gfpMatrix &res)
{
    // O(n^2)
    RandomShare::generate_active_open_randoms_online_offline(res);
}

void MultVerifier::random(gfpMatrix &res)
{
    ShareBundle r(res.rows(), res.cols());
    r.random();
    res = r.shares;
}

gfpMatrix MultVerifier::transform(vector<gfpScalar> input)
{
    gfpMatrix out(input.size(), 1);
    for (int i = 0; i < input.size(); i++) {
        out(i) = input.at(i);
    }
    return out;
}

gfpScalar MultVerifier::dot_product(gfpMatrix left, gfpMatrix right)
{
    assert(left.size() == right.size());
    gfpMatrix res = left.transpose() * right;
    return res(0);
}

// Batch dot_product: in left and right, each column is a vector.
void MultVerifier::dot_product_opt(const gfpMatrix &left, const gfpMatrix &right, gfpMatrix &res)
{
    assert(left.rows() == right.rows());
    assert(left.cols() == right.cols());
    assert(left.cols() == res.rows());
    res = (left.array() * right.array()).colwise().sum().transpose();
}

void MultVerifier::dot_product_opt(gfpMatrix &input, size_t colNum)
{
    size_t l = (input.rows() - 1) / 2;
    input(2*l, seqN(0, colNum)) = (input(seqN(0, l), seqN(0, colNum)).array() * input(seqN(l, l), seqN(0, colNum)).array()).colwise().sum().transpose();
}

void MultVerifier::init_coefficients(gfpVector &res, gfpScalar base_value, size_t pow_index)
{
    size_t len = res.size();
    for(size_t i = 0; i < len; i++) {
        res(i) = base_value.pow(pow_index);
        pow_index++;
    }
}


void MultVerifier::print_triples_num()
{
    if(offline_flag) {
        cout<<"The number of triples in offline phase is: "<<triples_num_offline<<endl;
    } else {
        cout<<"The number of triples in online phase is: "<<triples_num<<endl;
    }
}

// Push some random correct triples for debug.
void MultVerifier::test_verifier()
{
}

}