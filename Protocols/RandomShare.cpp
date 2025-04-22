#include "Protocols/RandomShare.h"
#include "Protocols/Bit.h"
#include "Protocols/MultVerifier.h"
#include "Protocols/InputVerifier.h"
#include <vector>

using Eigen::seq;
using Eigen::seqN;
using Eigen::last;
namespace hmmpc
{

// The queue store the preprocessed random share, which is nothing to do with the input of the party.
queue<gfpScalar> RandomShare::queueRandom; // [r]_t
queue<gfpScalar> RandomShare::queueRandomBit;

queue<gfpScalar> DoubleRandom::queueReducedRandom; // [r]_t, [r]_2t
queue<gfpScalar> DoubleRandom::queueTruncatedRandom; // [r/2^d]_t, [r]_t, [r_msb]_t
queue<gfpScalar> DoubleRandom::queueReducedTruncatedRandom; // [r/2^d]_t, [r]_2t, [r]_t, [r_msb]_t
queue<gfpScalar> DoubleRandom::queueUnboundedMultRandom; // ([b1], [b1^-1]), ([bi], [bi-1 * bi^-1]) for i = 2, ..., l. 
queue<gfpScalar> DoubleRandom::queueTruncatedRandomInML;
queue<gfpScalar> DoubleRandom::queueReducedTruncatedInML; // [r/2^d]_t, [r]_2t, [r]_t
queue<gfpScalar> DoubleRandom::queueReducedTruncatedWithPrecisionRandom; // [r/2^d]_t, [r]_2t, [r]_t
queue<gfpScalar> DoubleRandom::queueMultReducedRandomPair;
queue<gfpMatrix> DoubleRandom::queueReusedMultReducedRandomPair;
queue<gfpMatrix> DoubleRandom::queueReusedMultReducedRandomPair16;
queue<MatrixTriple> DoubleRandom::queueBeaverTriples;
/************************************************************************
 * 
 *       Definition of static member functions about RandomShare
 * 
 * **********************************************************************/
void RandomShare::get_random(queue<gfpScalar> &Q, gfpScalar &res)
{
    assert(Q.size()>0);
    res = Q.front();
    Q.pop();
    return;
}

void RandomShare::get_randoms(queue<gfpScalar> &Q, gfpMatrix &res)
{
    size_t num = res.rows()*res.cols();
    assert(num <= Q.size());
    for(size_t i = 0; i < num; i++){
        res(i) = Q.front();
        Q.pop();
    }
    return;
}

/**
 * @brief Generate num (at least num) random sharings using DN07 RandomShare Protocol.
 * Every execution of the RandomShare Protocol will generate (n-t) random sharings where n is number of players and t is the threshold.
 * 
 * @param num 
 */
void RandomShare::generate_random_sharings(size_t num)
{
    size_t bundle_size = bundles(num); // ceil( num / (n_player - threshold ))

    ShareBundle crude_inputs(n_players, bundle_size);// Each row corresponds the input of one party.
    // Input random sharing respectively from each party.
    vector<ShareBundle> individual_input(n_players);
    octetStreams os_send(n_players), os_receive(n_players);
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible. 
        individual_input[i].input_from_random_request(i, os_send, os_receive[i]);
    }

    // The multi-thread version.
    // The distributer should wait to send all the sharings.
    // The other receivers should wait to receive the shares.
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].finish_input_from(i, os_send, os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }
    
    // Use vandermonde matrix to extract randomness.
    gfpMatrix extracted_random = vandermonde_n_t.transpose()*crude_inputs.shares;
    for(size_t i = 0; i < extracted_random.size(); i++){
        queueRandom.push(extracted_random(i));
    }
}

void RandomShare::generate_random_sharings_active(size_t num)
{
    gfpMatrix passive_random_sharings(num+1, 1);
    generate_passive_random_sharings(passive_random_sharings);

    // Correctness check
    gfpMatrix active_open_randoms(num, 1);  // TODO_ZYS: seed+PRF
    generate_active_open_randoms(active_open_randoms);
    gfpMatrix coefficients(num+1, 1);
    // coefficients.block(0, 0, num, 1) = active_open_randoms;
    coefficients(seqN(0, num), 0) = active_open_randoms;
    coefficients(num) = gfpScalar(1);
    gfpScalar v = (coefficients.array() * passive_random_sharings.array()).sum();
    Share vShare;
    vShare.set_share(v);
    vShare.reveal_with_check();  // O(n^2)
    for(size_t i = 0; i < num; i++){
        queueRandom.push(passive_random_sharings(i));
    }
}


/**
 * @brief Generate a bundle of sharings of bits into the queueRandomBit queue.
 * 
 * @param num 
 */
void RandomShare::generate_random_bits(size_t num)
{
    ShareBundle R(num, 1);
    R.random();
    gfpMatrix &r = R.shares;
    
    ShareBundle r_square(num, 1);
    r_square.shares = square(r.array()); // *BUG LOG: We cannot use r.array().square() since we still use the origin r.
    
    // TODO: Modify with peusodo-random sharing of zero. 
    // In temporary, we use this unsecure operation but with the same complexity. (since PRSZ costs free)
    r_square.double_degree();
    gfpMatrix s = r_square.reveal(true);

    // r_square.reduce_degree(); // * BUG LOG: Remember to reduce degree when multiplying two shares of degree t. 
    // gfpMatrix s = r_square.reveal();
    for(size_t i = 0; i < s.size(); i++){
        if(s(i) == 0){ // TODO resample
            throw abort_error("[RandomShare]: get_random_bits samples 0");
        }
    }
    
    // gfpMatrix r_prime = s.array().rsqrt();
    // Use Batch Inversion in rsqrt
    gfpMatrix r_sqrt = s.array().sqrt();
    gfpMatrix r_prime(r_sqrt.rows(), r_sqrt.cols());
    batch_inversion(r_sqrt, r_prime);

    gfpMatrix res = ((r.array() * r_prime.array() ) + 1)/2;
    for(size_t i = 0; i < num; i++){
        queueRandomBit.push(res(i));
    }
    
    return;
}

void RandomShare::generate_random_bits_active(size_t num)
{
    ShareBundle R(num, 1);
    R.random();
    gfpMatrix &r = R.shares;
    
    ShareBundle r_square(num, 1);
    r_square.shares = square(r.array()); // *BUG LOG: We cannot use r.array().square() since we still use the origin r.
    
    // TODO: Modify with peusodo-random sharing of zero. 
    // In temporary, we use this unsecure operation but with the same complexity. (since PRSZ costs free)
    r_square.double_degree();
    gfpMatrix s = r_square.reveal(true);

    MultVerifier::push_triples(r, r, s, true);

    // r_square.reduce_degree(); // * BUG LOG: Remember to reduce degree when multiplying two shares of degree t. 
    // gfpMatrix s = r_square.reveal();
    for(size_t i = 0; i < s.size(); i++){
        if(s(i) == 0){ // TODO resample
            throw abort_error("[RandomShare]: get_random_bits samples 0");
        }
    }
    
    // gfpMatrix r_prime = s.array().rsqrt();
    // Use Batch Inversion in rsqrt
    gfpMatrix r_sqrt = s.array().sqrt();
    gfpMatrix r_prime(r_sqrt.rows(), r_sqrt.cols());
    batch_inversion(r_sqrt, r_prime);

    gfpMatrix res = ((r.array() * r_prime.array() ) + 1)/2;
    for(size_t i = 0; i < num; i++){
        queueRandomBit.push(res(i));
    }
    
    return;
}

// Replace with input_from_random_
void RandomShare::generate_random_sharings_PRG(size_t num)
{
    size_t bundle_size = bundles(num); // ceil( num / (n_player - threshold ))

    ShareBundle crude_inputs(n_players, bundle_size);// Each row corresponds the input of one party.
    
    // Input random sharing respectively from each party.
    vector<ShareBundle> individual_input(n_players);
    octetStreams os_send(n_party_PRG()), os_receive(n_players);
    vector<gfpMatrix> shares_prng(n_players);

    for(int i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible. 
        individual_input[i].input_from_random_request_PRG(i, os_send, shares_prng[i], os_receive[i]);
    }

    // The multi-thread version.
    // The distributer should wait to send all the sharings.
    // The other receivers should wait to receive the shares.
    for(int i = 0; i < n_players; i++){
        individual_input[i].finish_input_from_PRG(i, os_send, shares_prng[i], os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }

    // Use vandermonde matrix to extract randomness.
    gfpMatrix extracted_random = vandermonde_n_t.transpose()*crude_inputs.shares;
    for(size_t i = 0; i < extracted_random.size(); i++){
        queueRandom.push(extracted_random(i));
    }
    return;
}

void RandomShare::generate_random_sharings_PRG_active(size_t num)
{
    gfpMatrix passive_random_sharings(num, 1);
    generate_passive_random_sharings_PRG(passive_random_sharings);

    InputVerifier::push_matrix(passive_random_sharings, true);
    for(size_t i = 0; i < num; i++){
        queueRandom.push(passive_random_sharings(i));
    }
    return;
}

void RandomShare::generate_passive_random_sharings(gfpMatrix &res)
{
    size_t num = res.size();
    size_t bundle_size = bundles(num); // ceil( num / (n_player - threshold ))

    ShareBundle crude_inputs(n_players, bundle_size);// Each row corresponds the input of one party.
    // Input random sharing respectively from each party.
    vector<ShareBundle> individual_input(n_players);
    octetStreams os_send(n_players), os_receive(n_players);
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible. 
        individual_input[i].input_from_random_request(i, os_send, os_receive[i]);
    }

    // The multi-thread version.
    // The distributer should wait to send all the sharings.
    // The other receivers should wait to receive the shares.
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].finish_input_from(i, os_send, os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }
    
    // Use vandermonde matrix to extract randomness.
    // DEBUG LOG: Eigen::Matrix.resize() could change the value in the matrix.
    gfpMatrix extracted_random = vandermonde_n_t.transpose()*crude_inputs.shares;
    for (int i = 0; i < res.size(); i++)
    {
        res(i) = extracted_random(i);
    }
}

void RandomShare::generate_active_random_sharings(gfpMatrix &res)
{
    gfpMatrix passive_random_sharings(2*res.size()+1, 1);
    generate_passive_random_sharings(passive_random_sharings);

    // generate active open randoms
    ShareBundle active_open_randoms(res.size(), 1);
    active_open_randoms.shares = passive_random_sharings(seqN(res.size()+1, res.size()), seqN(0, 1));
    active_open_randoms.reveal_with_check();

    // Correctness check
    gfpMatrix coefficients(res.size()+1, 1);
    // coefficients.block(0, 0, res.size(), 1) = active_open_randoms;
    coefficients(seqN(0, res.size()), 0) = active_open_randoms.get_secrets();
    coefficients(res.size()) = gfpScalar(1);
    Share vShare;
    vShare.share = (coefficients.array() * passive_random_sharings(seqN(0, res.size()+1), 0).array()).sum();
    vShare.reveal_with_check();  // O(n^2)
    res = passive_random_sharings(seqN(0, res.size()), 0).reshaped<Eigen::RowMajor>(res.rows(), res.cols()); 
}

void RandomShare::generate_active_open_randoms(gfpMatrix &res)
{
    // Can not be equal to zero in paper [LN17], but zero also can satisfy the statistical security.
    ShareBundle random_data(res.rows(), res.cols());
    generate_passive_random_sharings(random_data.shares);
    res = random_data.reveal_with_check();  // O(n^2)
}

void RandomShare::generate_passive_random_sharings_PRG(gfpMatrix &res)
{
    size_t num = res.size();
    size_t bundle_size = bundles(num); // ceil( num / (n_player - threshold ))

    ShareBundle crude_inputs(n_players, bundle_size);// Each row corresponds the input of one party.
    
    // Input random sharing respectively from each party.
    vector<ShareBundle> individual_input(n_players);
    octetStreams os_send(n_party_PRG()), os_receive(n_players);
    vector<gfpMatrix> shares_prng(n_players);

    for(int i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible. 
        individual_input[i].input_from_random_request_PRG(i, os_send, shares_prng[i], os_receive[i]);
    }

    // The multi-thread version.
    // The distributer should wait to send all the sharings.
    // The other receivers should wait to receive the shares.
    for(int i = 0; i < n_players; i++){
        individual_input[i].finish_input_from_PRG(i, os_send, shares_prng[i], os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }

    // Use vandermonde matrix to extract randomness.
    gfpMatrix extracted_random = vandermonde_n_t.transpose()*crude_inputs.shares;
    for (int i = 0; i < res.size(); i++)
    {
        res(i) = extracted_random(i);
    }
}

void RandomShare::generate_active_random_sharings_PRG(gfpMatrix &res)
{
    gfpMatrix passive_random_sharings(2*res.size()+1, 1);
    generate_passive_random_sharings_PRG(passive_random_sharings);

    // generate active open randoms
    ShareBundle active_open_randoms(res.size(), 1);
    active_open_randoms.shares = passive_random_sharings(seqN(res.size()+1, res.size()), seqN(0, 1));
    active_open_randoms.reveal_with_check();

    // Correctness check
    gfpMatrix coefficients(res.size()+1, 1);
    // coefficients.block(0, 0, res.size(), 1) = active_open_randoms;
    coefficients(seqN(0, res.size()), 0) = active_open_randoms.get_secrets();
    coefficients(res.size()) = gfpScalar(1);
    Share vShare;
    vShare.share = (coefficients.array() * passive_random_sharings(seqN(0, res.size()+1), 0).array()).sum();
    vShare.reveal_with_check();  // O(n^2)
    res = passive_random_sharings(seqN(0, res.size()), 0).reshaped<Eigen::RowMajor>(res.rows(), res.cols()); 
}

void RandomShare::generate_active_open_randoms_PRG(gfpMatrix &res)
{
    // Can not be equal to zero in paper [LN17], but zero also can satisfy the statistical security.
    ShareBundle random_data(res.rows(), res.cols());
    generate_passive_random_sharings_PRG(random_data.shares);
    res = random_data.reveal_with_check();  // O(n^2)
}

void RandomShare::generate_active_open_randoms_online_offline(gfpMatrix &res)
{
    ShareBundle random_data(res.rows(), res.cols());
    random_data.random();  // Note_ZYS: execute in the offline phase which could provide active security. But passive security is OK.
    res = random_data.reveal_with_check();  // O(n^2)
}

/************************************************************************
 * 
 *       Definition of static member functions about DoubleRandom
 * 
 * **********************************************************************/

void DoubleRandom::get_random_pair(queue<gfpScalar>&Q, gfpScalar &r, gfpScalar &aux_r)
{
    assert(Q.size()>=2);
    r = Q.front();
    Q.pop();
    aux_r = Q.front();
    Q.pop();
    return;
}


void DoubleRandom::get_random_pairs(queue<gfpScalar>&Q, gfpMatrix &r, gfpMatrix &aux_r)
{
    assert(r.size()*2 <= Q.size());
    for(size_t i = 0; i < r.size(); i++){
        r(i) = Q.front();
        Q.pop();
        aux_r(i) = Q.front();
        Q.pop();
    }
    return;
}

void DoubleRandom::get_random_triple(queue<gfpScalar>&Q, gfpScalar &r, gfpScalar &aux_r, gfpScalar &sub_r)
{
    assert(Q.size()>=3);
    r = Q.front();
    Q.pop();
    aux_r = Q.front();
    Q.pop();
    sub_r = Q.front();
    Q.pop();
    return;
}

void DoubleRandom::get_random_triples(queue<gfpScalar> &Q, gfpMatrix &r, gfpMatrix &aux_r, gfpMatrix &sub_r)
{
    assert(r.size()*3 <= Q.size());
    for(size_t i = 0; i < r.size(); i++){
        r(i) = Q.front();
        Q.pop();
        aux_r(i) = Q.front();
        Q.pop();
        sub_r(i) = Q.front();
        Q.pop();
    }
    return;
}

void DoubleRandom::get_random_quadruple(queue<gfpScalar>&Q, gfpScalar &r, gfpScalar &aux_r, gfpScalar &aux_aux_r, gfpScalar &sub_r)
{
    assert(Q.size()>=4);
    r = Q.front();
    Q.pop();
    aux_r = Q.front();
    Q.pop();
    aux_aux_r = Q.front();
    Q.pop();
    sub_r = Q.front();
    Q.pop();
    return;
}

void DoubleRandom::get_random_quadruples(queue<gfpScalar> &Q, gfpMatrix &r, gfpMatrix &aux_r, gfpMatrix &aux_aux_r, gfpMatrix &sub_r)
{
    assert(r.size()*4 <= Q.size());
    for(size_t i = 0; i < r.size(); i++){
        r(i) = Q.front();
        Q.pop();
        aux_r(i) = Q.front();
        Q.pop();
        aux_aux_r(i) = Q.front();
        Q.pop();
        sub_r(i) = Q.front();
        Q.pop();
    }
    return;
}

void DoubleRandom::get_mult_reduced_random_sharing_pairs(gfpMatrix &r1, gfpMatrix &aux_r1, gfpMatrix &r2, gfpMatrix &aux_r2, gfpMatrix &prod)
{
    assert(r1.size()*5 <= queueMultReducedRandomPair.size());
    for(size_t i = 0; i < r1.size(); i++){
        r1(i) = queueMultReducedRandomPair.front();
        queueMultReducedRandomPair.pop();
        aux_r1(i) = queueMultReducedRandomPair.front();
        queueMultReducedRandomPair.pop();
        r2(i) = queueMultReducedRandomPair.front();
        queueMultReducedRandomPair.pop();
        aux_r2(i) = queueMultReducedRandomPair.front();
        queueMultReducedRandomPair.pop();
        prod(i) = queueMultReducedRandomPair.front();
        queueMultReducedRandomPair.pop();
    }
    return;
}

void DoubleRandom::get_reused_mult_reduced_random_sharing_pairs(gfpMatrix &r1, gfpMatrix &aux_r1, gfpMatrix &r2, gfpMatrix &aux_r2, gfpMatrix &prod)
{
    assert(r1.size() <= queueReusedMultReducedRandomPair.size());
    for(size_t i = 0; i < r1.rows(); i++){
        for(size_t j = 0; j < r1.cols(); j++) {
            gfpMatrix tmp = queueReusedMultReducedRandomPair.front();
            queueReusedMultReducedRandomPair.pop();
            // ([r]_t, [r]_2t)
            r1(i, j) = tmp(0, 0);
            aux_r1(i, j) = tmp(1, 0);
            for (size_t k = 0; k < 4; k++)
            {
                // ([r1]_t, [r1]_2t),...,([r4]_t, [r4]_2t)
                size_t index = k * r1.rows();
                r2(i+index, j) = tmp(0, k+1);
                aux_r2(i+index, j) = tmp(1, k+1);
            }
            // [r*r1]_t
            prod(i, j) = tmp(0, 5);
            // [r*r2]_t
            prod(i+r1.rows(), j) = tmp(1, 5);
            // [r*r3]_t
            prod(i+2*r1.rows(), j) = tmp(0, 6);
            // [r*r4]_t
            prod(i+3*r1.rows(), j) = tmp(1, 6);
        }
    }
    return;
}

void DoubleRandom::get_reused_mult_reduced_random_sharing_pairs_16(gfpMatrix &r1, gfpMatrix &aux_r1, gfpMatrix &r2, gfpMatrix &aux_r2, gfpMatrix &prod)
{
    assert(r1.size() <= queueReusedMultReducedRandomPair16.size());
    for(size_t i = 0; i < r1.rows(); i++){
        for(size_t j = 0; j < r1.cols(); j++) {
            gfpMatrix tmp = queueReusedMultReducedRandomPair16.front();
            queueReusedMultReducedRandomPair16.pop();
            // ([r]_t, [r]_2t)
            r1(i, j) = tmp(0, 0);
            aux_r1(i, j) = tmp(1, 0);
            for (size_t k = 0; k < 16; k++)
            {
                // ([r1]_t, [r1]_2t),...,([r16]_t, [r16]_2t)
                size_t index = k * r1.rows();
                r2(i+index, j) = tmp(0, k+1);
                aux_r2(i+index, j) = tmp(1, k+1);
            }
            for (size_t k = 0, m = 0; k < 16; k++, m++)
            {
                // [r*r1]_t,...,[r*r16]_t
                prod(i+(k++)*r1.rows(), j) = tmp(0, 17+m);
                prod(i+k*r1.rows(), j) = tmp(1, 17+m);
            }
        }
    }
    return;
}

void DoubleRandom::get_beaver_triple(gfpMatrix &a, gfpMatrix &b, gfpMatrix &c)
{
    assert(queueBeaverTriples.size() > 0);
    MatrixTriple mat = queueBeaverTriples.front();
    assert(a.rows() == mat.left.rows());
    assert(a.cols() == mat.left.cols());
    assert(b.rows() == mat.right.rows());
    assert(b.cols() == mat.right.cols());
    assert(c.rows() == mat.product.rows());
    assert(c.cols() == mat.product.cols());

    a = mat.left;
    b = mat.right;
    c = mat.product;
    queueBeaverTriples.pop();

    return;
}

/**********************************************************
 * *      Reduced Random Sharings - ( [r]_t, [r]_2t )
 * *********************************************************/
/**
 * @brief Generate num (at least num) double random sharings using DN07 DoubleRandom Protocol.
 * Every execution of the DoubleRandom Protocol will generate (n-t) double random sharings where n is number of players and t is the threshold.
 * Each pair of double random sharing consists of ([r]_t, [r]_2t)
 * 
 * @param num 
 */
void DoubleRandom::generate_reduced_random_sharings(size_t num)
{
    size_t bundle_size = bundles(num);
    DoubleShareBundle crude_inputs(n_players, bundle_size); // Each row corresponds to one party's input
    // Input reduced random sharing respectively from each party
    vector<DoubleShareBundle> individual_input(n_players);
    octetStreams os_send(n_players), os_receive(n_players);
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible. 
        individual_input[i].input_from_random_request(i, os_send, os_receive[i]);
    }

    // The multi-thread version.
    // The distributer should wait to send all the sharings.
    // The other receivers should wait to receive the shares.
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].finish_input_from(i, os_send, os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }

    // We can calculate extracted_t_random and extracted_2t_random in one matrix multiplication.
    // *Original matrix multiplication
    // gfpMatrix extracted_t_random = vandermonde_n_t.transpose() * crude_inputs.shares;
    // gfpMatrix extracted_2t_random = vandermonde_n_t.transpose() * crude_inputs.aux_shares; 
    // for(size_t i = 0; i < extracted_t_random.size(); i++){
    //     queueReducedRandom.push(extracted_t_random(i));
    //     queueReducedRandom.push(extracted_2t_random(i));
    // }


    // Use vandermonde matrix to extract randomness
    // We incorporate two matrix multiplication into one matrix multiplication.
    gfpMatrix extracted_random = vandermonde_n_t.transpose() * 
        (gfpMatrix(n_players, bundle_size<<1)<<crude_inputs.shares, crude_inputs.aux_shares).finished();
    
    for(size_t i = 0; i < extracted_random.rows(); i++){
        // j indicates the t-sharing and k indicates the 2t-sharing
        for(size_t j = 0, k = bundle_size; k < extracted_random.cols(); j++, k++){
            queueReducedRandom.push(extracted_random(i, j));
            queueReducedRandom.push(extracted_random(i, k));
        }
    }
    
}

void DoubleRandom::generate_reduced_random_sharings_active(size_t num)
{
    size_t bundle_size = bundles(num);
    DoubleShareBundle crude_inputs(n_players, bundle_size); // Each row corresponds to one party's input
    // Input reduced random sharing respectively from each party
    vector<DoubleShareBundle> individual_input(n_players);
    octetStreams os_send(n_players), os_receive(n_players);
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible. 
        individual_input[i].input_from_random_request(i, os_send, os_receive[i]);
    }

    // The multi-thread version.
    // The distributer should wait to send all the sharings.
    // The other receivers should wait to receive the shares.
    for(size_t i = 0; i < n_players; i++){
        individual_input[i].finish_input_from(i, os_send, os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }

    // We can calculate extracted_t_random and extracted_2t_random in one matrix multiplication.
    // *Original matrix multiplication
    // gfpMatrix extracted_t_random = vandermonde_n_t.transpose() * crude_inputs.shares;
    // gfpMatrix extracted_2t_random = vandermonde_n_t.transpose() * crude_inputs.aux_shares; 
    // for(size_t i = 0; i < extracted_t_random.size(); i++){
    //     queueReducedRandom.push(extracted_t_random(i));
    //     queueReducedRandom.push(extracted_2t_random(i));
    // }


    // Use vandermonde matrix to extract randomness
    // We incorporate two matrix multiplication into one matrix multiplication.
    gfpMatrix extracted_random = vandermonde_n_t.transpose() * 
        (gfpMatrix(n_players, bundle_size<<1)<<crude_inputs.shares, crude_inputs.aux_shares).finished();

    // Check the correctness of t-sharing randoms
    // gfpMatrix padding_checked = extracted_random.block(0, 0, extracted_random.rows(), bundle_size).reshaped<Eigen::RowMajor>();
    gfpMatrix padding_checked = extracted_random(seqN(0, extracted_random.rows()), seqN(0, bundle_size)).reshaped<Eigen::RowMajor>();
    gfpMatrix active_open_randoms(padding_checked.size(), 1);  // TODO: seed+PRF
    generate_active_open_randoms(active_open_randoms);
    gfpMatrix r(1, 1);
    generate_passive_random_sharings(r);
    Share vShare;
    vShare.set_share((active_open_randoms.array()*padding_checked.array()).sum() + r(0, 0));
    vShare.reveal_with_check();  // O(n^2)
    
    for(size_t i = 0; i < extracted_random.rows(); i++){
        // j indicates the t-sharing and k indicates the 2t-sharing
        for(size_t j = 0, k = bundle_size; k < extracted_random.cols(); j++, k++){
            queueReducedRandom.push(extracted_random(i, j));
            queueReducedRandom.push(extracted_random(i, k));
        }
    }
}

// Replace with input_from_random_request_PRG.
void DoubleRandom::generate_reduced_random_sharings_PRG(size_t num)
{
    size_t bundle_size = bundles(num);
    DoubleShareBundle crude_inputs(n_players, bundle_size); // Each row corresponds to one party's input
    // Input reduced random sharing respectively from each party

    // * The part that is replaced.
    vector<DoubleShareBundle> individual_input(n_players);
    octetStreams os_send(n_players), os_receive(n_players);
    vector<gfpMatrix>shares_prng(n_players);

    for(int i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible.
        individual_input[i].input_from_random_request_PRG(i, os_send, shares_prng[i], os_receive[i]); 
    }

    for(int i = 0; i < n_players; i++){
        individual_input[i].input_aux_from_random_PRG(i);
    }

    for(int i = 0; i < n_players; i++){
        individual_input[i].finish_input_from_PRG(i, os_send, shares_prng[i], os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }
    // *

    gfpMatrix extracted_random = vandermonde_n_t.transpose() * 
        (gfpMatrix(n_players, bundle_size<<1)<<crude_inputs.shares, crude_inputs.aux_shares).finished();
    
    for(size_t i = 0; i < extracted_random.rows(); i++){
        // j indicates the t-sharing and k indicates the 2t-sharing
        for(size_t j = 0, k = bundle_size; k < extracted_random.cols(); j++, k++){
            queueReducedRandom.push(extracted_random(i, j));
            queueReducedRandom.push(extracted_random(i, k));
        }
    }
}

void DoubleRandom::generate_reduced_random_sharings_PRG_active(size_t num)
{
    size_t bundle_size = bundles(num);
    DoubleShareBundle crude_inputs(n_players, bundle_size); // Each row corresponds to one party's input
    // Input reduced random sharing respectively from each party

    // * The part that is replaced.
    vector<DoubleShareBundle> individual_input(n_players);
    octetStreams os_send(n_players), os_receive(n_players);
    vector<gfpMatrix>shares_prng(n_players);

    for(int i = 0; i < n_players; i++){
        individual_input[i].resize(1, bundle_size);
        // Each party only requests to distribute the sharings 
        // so that all the parties could start to distribute sharings as soon as possible.
        individual_input[i].input_from_random_request_PRG(i, os_send, shares_prng[i], os_receive[i]); 
    }

    for(int i = 0; i < n_players; i++){
        individual_input[i].input_aux_from_random_PRG(i);
    }

    for(int i = 0; i < n_players; i++){
        individual_input[i].finish_input_from_PRG(i, os_send, shares_prng[i], os_receive[i]);
        crude_inputs.set_row(i, individual_input[i]);
    }
    // *

    gfpMatrix extracted_random = vandermonde_n_t.transpose() * 
        (gfpMatrix(n_players, bundle_size<<1)<<crude_inputs.shares, crude_inputs.aux_shares).finished();

    // Check the correctness of t-sharing randoms
    // gfpMatrix padding_checked = extracted_random.block(0, 0, extracted_random.rows(), bundle_size).reshaped<Eigen::RowMajor>();
    gfpMatrix padding_checked = extracted_random(seqN(0, extracted_random.rows()), seqN(0, bundle_size)).reshaped<Eigen::RowMajor>();
    InputVerifier::push_matrix(padding_checked, true);
    
    for(size_t i = 0; i < extracted_random.rows(); i++){
        // j indicates the t-sharing and k indicates the 2t-sharing
        for(size_t j = 0, k = bundle_size; k < extracted_random.cols(); j++, k++){
            queueReducedRandom.push(extracted_random(i, j));
            queueReducedRandom.push(extracted_random(i, k));
        }
    }
}

/**********************************************************
 * *      Truncated Random Sharings - ( [r/2^d]_t, [r]_t , [r_{msb}]_t)
 * *********************************************************/
void DoubleRandom::generate_truncated_random_sharings(size_t num)
{
    BitBundle bits(num, BITS_LENGTH);
    bits.random();
    
    gfpMatrix trunc_bits(num, BITS_LENGTH);
    trunc_bits.leftCols(INT_PRECISION) = bits.shares.rightCols(INT_PRECISION);
    for(size_t i = INT_PRECISION; i < BITS_LENGTH; i++){
        // Fill the empty bits with the MSB of the original bits
        trunc_bits.col(i) = bits.shares.col(BITS_LENGTH - 1);
    }

    // *Original matrix multiplication
    // gfpMatrix res = bits.shares * bits_coeff;
    // gfpMatrix res_trunc = trunc_bits * bits_coeff; 

    // We incorporate two matrix multiplication into one matrix multiplication.
    gfpMatrix res = (gfpMatrix(num<<1, BITS_LENGTH)<< trunc_bits, bits.shares).finished() * bits_coeff;
    for(size_t i = 0, j = num; i < num; i++, j++){
        queueTruncatedRandom.push(res(i));
        queueTruncatedRandom.push(res(j));
        queueTruncatedRandom.push(bits.shares(i, BITS_LENGTH-1));// msb
    }
}

// Truncate with specific precision
void DoubleRandom::generate_truncated_random_sharings(size_t num, size_t precision)
{
    BitBundle bits(num, BITS_LENGTH);
    bits.random();
    
    size_t int_precision = BITS_LENGTH - precision;
    gfpMatrix trunc_bits(num, BITS_LENGTH);
    trunc_bits.leftCols(int_precision) = bits.shares.rightCols(int_precision);
    for(size_t i = int_precision; i < BITS_LENGTH; i++){
        // Fill the empty bits with the MSB of the original bits
        trunc_bits.col(i) = bits.shares.col(BITS_LENGTH - 1);
    }

    // We incorporate two matrix multiplication into one matrix multiplication.
    gfpMatrix res = (gfpMatrix(num<<1, BITS_LENGTH)<< trunc_bits, bits.shares).finished() * bits_coeff;
    for(size_t i = 0, j = num; i < num; i++, j++){
        queueTruncatedRandomInML.push(res(i));
        queueTruncatedRandomInML.push(res(j));
    }
}


/**********************************************************************
 * *      Reduced Truncated Random Sharings - ( [r/2^d]_t, [r]_2t, [r]_t, [r_msb]_t)
 * *********************************************************************/
void DoubleRandom::generate_reduced_truncated_random_sharings(size_t num)
{
    BitBundle bitsBundle(num, BITS_LENGTH);
    bitsBundle.random();
    gfpMatrix bits = bitsBundle.shares;
    gfpMatrix trunc_bits(num, BITS_LENGTH);
    trunc_bits.leftCols(INT_PRECISION) = bits.rightCols(INT_PRECISION);
    for(size_t i = INT_PRECISION; i < BITS_LENGTH; i++){
        // Fill the empty bits with the MSB of the original bits
        trunc_bits.col(i) = bits.col(BITS_LENGTH - 1);
    }
    bits = bits.array().square(); // Each bit turns to 2t-sharing
    gfpMatrix res = (gfpMatrix(num*3, BITS_LENGTH)<< trunc_bits, bits, bitsBundle.shares).finished() * bits_coeff;

    for(size_t i = 0, j = num, k = 2*num; i < num; i++, j++, k++){
        queueReducedTruncatedRandom.push(res(i));
        queueReducedTruncatedRandom.push(res(j));
        queueReducedTruncatedRandom.push(res(k));
        queueReducedTruncatedRandom.push(trunc_bits(i, BITS_LENGTH-1));// msb
    }
}

void DoubleRandom::generate_reduced_truncated_random_sharings(queue<gfpScalar> &Q, size_t num, size_t precision)
{
    BitBundle R(num, BITS_LENGTH);
    R.random();
    gfpMatrix bits = R.shares;

    gfpMatrix trunc_bits(num, BITS_LENGTH);
    size_t int_precision = BITS_LENGTH - precision;
    trunc_bits.leftCols(int_precision) = bits.rightCols(int_precision);
    for(size_t i = int_precision; i < BITS_LENGTH; i++){
        // Fill the empty bits with the MSB of the original bits
        trunc_bits.col(i) = bits.col(BITS_LENGTH - 1);
    }

    bits = bits.array().square(); // Each bit turns to 2t-sharing
    gfpMatrix res = (gfpMatrix(num*3, BITS_LENGTH)<< trunc_bits, bits, R.shares).finished() * bits_coeff;

    for(size_t i = 0, j = num, k = 2*num; i < num; i++, j++, k++){
        Q.push(res(i));
        Q.push(res(j));
        Q.push(res(k));
    }
}

// Generate reduced and truncated random sharings with different precision.
// num_repetitions: the number of repetitions
void DoubleRandom::generate_reduced_truncated_random_sharings(queue<gfpScalar> &Q, size_t num, vector<size_t> &precision, size_t num_repetitions)
{
    size_t total = num*num_repetitions;
    BitBundle R(total, BITS_LENGTH);
    R.random();
    gfpMatrix bits = R.shares;

    gfpMatrix trunc_bits(total, BITS_LENGTH);
    vector<size_t> int_precision;
    for(size_t i = 0; i < num_repetitions; i++){
        for(size_t j = 0; j < num; j++){
            int_precision.push_back(BITS_LENGTH - precision[j]);
        }
    }

    for(size_t i = 0; i < total; i++){
        trunc_bits.row(i).leftCols(int_precision[i]) = bits.row(i).rightCols(int_precision[i]);
        trunc_bits.row(i).rightCols(BITS_LENGTH - int_precision[i]).setConstant(bits(i, BITS_LENGTH-1));
    }

    bits = bits.array().square(); // Each bit turns to 2t-sharing
    gfpMatrix res = (gfpMatrix(total*3, BITS_LENGTH)<< trunc_bits, bits, R.shares).finished() * bits_coeff;

    for(size_t i = 0, j = total, k = 2*total; i < total; i++, j++, k++){
        Q.push(res(i));
        Q.push(res(j));
        Q.push(res(k));
    }
}


void DoubleRandom::generate_non_zero_random_column_wise(gfpMatrix &res, gfpMatrix &prod_reveal, bool active)
{
    assert(res.cols() == 2);
    ShareBundle R(res.rows(), 2);
    R.random();
    ShareBundle prod(res.rows(), 1);
    prod.shares = R.shares.col(0).array() * R.shares.col(1).array();
    prod.double_degree();
    prod.reveal(active);
    if(active) {
        // Verify multiplication & public
        MultVerifier::push_triples(R.shares.col(0), R.shares.col(1), prod.secret(), true);
    }
    res = R.shares;
    prod_reveal = prod.secret();
    vector<pair<int, int>> zero_index;
    size_t zero_num = find_zero(prod_reveal, zero_index);
    if(zero_num > 0) {
        gfpMatrix R_extend(zero_num, 2);
        gfpMatrix prod_extend(zero_num, 1);
        generate_non_zero_random_column_wise(R_extend, prod_extend, active);
        for(int i = 0; i < zero_num; i++) {
            res.row(zero_index[i].first) = R_extend.row(i);
            prod_reveal.row(zero_index[i].first) = prod_extend.row(i);
        }
    }
}

void DoubleRandom::generate_non_zero_random_row_wise(gfpMatrix &res, gfpMatrix &prod_reveal, bool active)
{
    int block_num = res.rows() / 2;
    assert(res.rows() % 2 == 0);
    ShareBundle R(res.rows(), res.cols());
    R.random();
    ShareBundle prod(block_num, res.cols());
    prod.shares = R.shares(seqN(0, block_num), seqN(0, res.cols())).array() * R.shares(seqN(block_num, block_num), seqN(0, res.cols())).array();
    prod.double_degree();
    prod.reveal(active);
    if(active) {
        // Verify multiplication & public
        MultVerifier::push_triples(R.shares(seqN(0, block_num), seqN(0, res.cols())), R.shares(seqN(block_num, block_num), seqN(0, res.cols())), prod.secret(), true);
    }
    res = R.shares;
    prod_reveal = prod.secret();
    vector<pair<int, int>> zero_index;
    size_t zero_num = find_zero(prod_reveal, zero_index);
    if(zero_num > 0) {
        gfpMatrix R_extend(zero_num, 2);
        gfpMatrix prod_extend(zero_num, 1);
        generate_non_zero_random_column_wise(R_extend, prod_extend, active);
        for(int i = 0; i < zero_num; i++) {
            res(zero_index[i].first, zero_index[i].second) = R_extend(i, 0);
            res(zero_index[i].first+block_num, zero_index[i].second) = R_extend(i, 1);
            prod_reveal(zero_index[i].first, zero_index[i].second) = prod_extend(i, 0);
        }
    }
}

/**********************************************************************
 * *      Random Sharings for Unbounded Multiplication 
 * An unbounded multiplication instance
 * ( [b1]_t , [b1^-1]_t )
 * ( [b2]_t , [b1 * b2^-1]_t,)
 * ( [b3]_t , [b2 * b3^-1]_t)
 * ...
 * ( [bi]_t , [bi-1 * bi^-1]_t) for i = 2, ..., l
 * *********************************************************************/
// One unbounded-mult instance.
// The following is unbounded multiplication in parallel, which is more practical.
void DoubleRandom::generate_unbounded_random_sharings(size_t num)
{
    // TODO_ZYS: R can not be zero since we want to compute the inversion.
    gfpMatrix rand_b(num, 2);
    gfpMatrix prod_reveal(num, 1);
    generate_non_zero_random_column_wise(rand_b, prod_reveal);
    gfpMatrix tmp_b_bPrime(num-1, 2); // (bi, bi')
    
    tmp_b_bPrime.col(0)(seqN(num, num-1)) = rand_b.col(0)(seqN(0, num-1));
    tmp_b_bPrime.col(1)(seqN(num, num-1)) = rand_b.col(1)(seqN(1, num-1));
    // tmp_b_bPrime(seqN(num, num-1), seq(0, 0)) = b_bPrime(seqN(0, num-1), seq(0, 0));
    // tmp_b_bPrime(seqN(num, num-1), seq(1, 1)) = b_bPrime(seqN(1, num-1), seq(1, 1));
    
    ShareBundle prod(num-1, 1); // Reduce the use of Double random.
    prod.shares = tmp_b_bPrime.col(0).array() * tmp_b_bPrime.col(1).array();
    prod.reduce_degree();

    // gfpMatrix B_inv = B.secret().array().inverse(); // B_inv 
    // *Use Batch Inversion
    gfpMatrix B_inv(prod_reveal.rows(), prod_reveal.cols());
    batch_inversion(prod_reveal, B_inv);

    gfpVector prod_b(num);
    prod_b(0) = rand_b(0, 1); // b0'
    prod_b(seqN(1, num-1)) = prod.shares;

    gfpVector b_bPrime = B_inv.array() * prod_b.array();
    for(size_t i = 0; i < num; i++){
        queueUnboundedMultRandom.push(rand_b(i, 0));
        queueUnboundedMultRandom.push(b_bPrime(i));
    }
}

void DoubleRandom::generate_unbounded_random_sharings_active(size_t num)
{
    // TODO_ZYS: R can not be zero since we want to compute the inversion.
    gfpMatrix rand_b(num, 2);
    gfpMatrix prod_reveal(num, 1);
    generate_non_zero_random_column_wise(rand_b, prod_reveal, true);
    gfpMatrix tmp_b_bPrime(num-1, 2); // (bi, bi')
    
    tmp_b_bPrime.col(0)(seqN(num, num-1)) = rand_b.col(0)(seqN(0, num-1));
    tmp_b_bPrime.col(1)(seqN(num, num-1)) = rand_b.col(1)(seqN(1, num-1));
    // tmp_b_bPrime(seqN(num, num-1), seq(0, 0)) = b_bPrime(seqN(0, num-1), seq(0, 0));
    // tmp_b_bPrime(seqN(num, num-1), seq(1, 1)) = b_bPrime(seqN(1, num-1), seq(1, 1));
    
    ShareBundle prod(num-1, 1); // Reduce the use of Double random.
    prod.shares = tmp_b_bPrime.col(0).array() * tmp_b_bPrime.col(1).array();
    prod.reduce_degree();
    // Verify
    MultVerifier::push_triples(tmp_b_bPrime.col(0), tmp_b_bPrime.col(1), prod.shares, true);

    // gfpMatrix B_inv = B.secret().array().inverse(); // B_inv 
    // *Use Batch Inversion
    gfpMatrix B_inv(prod_reveal.rows(), prod_reveal.cols());
    batch_inversion(prod_reveal, B_inv);

    gfpVector prod_b(num);
    prod_b(0) = rand_b(0, 1); // b0'
    prod_b(seqN(1, num-1)) = prod.shares;

    gfpVector b_bPrime = B_inv.array() * prod_b.array();
    for(size_t i = 0; i < num; i++){
        queueUnboundedMultRandom.push(rand_b(i, 0));
        queueUnboundedMultRandom.push(b_bPrime(i));
    }
}

// Each row corresponds to an unbounded multiplication instance.
void DoubleRandom::generate_unbounded_random_sharings(size_t xSize, size_t ySize)
{
    gfpMatrix rand(xSize<<1, ySize);
    gfpMatrix prod_reveal(xSize, ySize);
    generate_non_zero_random_row_wise(rand, prod_reveal);

    gfpMatrix rand_b(xSize, ySize-1), rand_b_prime(xSize, ySize-1);

    // rand_b = b | b
    // rand_b = b'| b'
    rand_b(seqN(0, xSize), seqN(0, ySize-1)) = rand(seqN(0, xSize), seqN(0, ySize-1));
    rand_b_prime(seqN(0, xSize), seqN(0, ySize-1)) = rand(seqN(xSize, xSize), seqN(1, ySize-1));

    ShareBundle prod(xSize, ySize-1);  // Reduce the use of Double random.
    prod.shares = rand_b.array() * rand_b_prime.array();
    prod.reduce_degree();

    // gfpMatrix B_inv =  B.secret().array().inverse(); // B_inv //TODO: Evaluate the ex_gcd and pow to implement inverse. 
                                                    // BUG LOG: inverse() dose not return the reference.
    // *Use Batch Inversion
    gfpMatrix B_inv(prod_reveal.rows(), prod_reveal.cols());
    batch_inversion(prod_reveal, B_inv);

    gfpMatrix X(xSize, ySize);
    X.col(0) = rand(seqN(xSize, xSize), 0); // b'
    X(seqN(0, xSize), seqN(1, ySize-1)) = prod.shares;
    
    gfpMatrix res(xSize, ySize);
    res = B_inv.array() * X.array();
    for(size_t i = 0; i < xSize; i++)
        for(size_t j = 0; j < ySize; j++){
            queueUnboundedMultRandom.push(rand(i, j));
            queueUnboundedMultRandom.push(res(i, j));
        }
    return;
}

void DoubleRandom::generate_unbounded_random_sharings_active(size_t xSize, size_t ySize)
{
    gfpMatrix rand(xSize<<1, ySize);
    gfpMatrix prod_reveal(xSize, ySize);
    generate_non_zero_random_row_wise(rand, prod_reveal, true);

    gfpMatrix rand_b(xSize, ySize-1), rand_b_prime(xSize, ySize-1);

    // rand_b = b | b
    // rand_b = b'| b'
    rand_b(seqN(0, xSize), seqN(0, ySize-1)) = rand(seqN(0, xSize), seqN(0, ySize-1));
    rand_b_prime(seqN(0, xSize), seqN(0, ySize-1)) = rand(seqN(xSize, xSize), seqN(1, ySize-1));

    ShareBundle prod(xSize, ySize-1); // Reduce the use of Double random.
    prod.shares = rand_b.array() * rand_b_prime.array();
    prod.reduce_degree();
    // Verify
    MultVerifier::push_triples(rand_b, rand_b_prime, prod.shares, true);

    // gfpMatrix B_inv =  B.secret().array().inverse(); // B_inv //TODO: Evaluate the ex_gcd and pow to implement inverse. 
                                                    // BUG LOG: inverse() dose not return the reference.
    // *Use Batch Inversion
    gfpMatrix B_inv(prod_reveal.rows(), prod_reveal.cols());
    batch_inversion(prod_reveal, B_inv);

    gfpMatrix X(xSize, ySize);
    X.col(0) = rand(seqN(xSize, xSize), 0); // b'
    X(seqN(0, xSize), seqN(1, ySize-1)) = prod.shares;
    
    gfpMatrix res(xSize, ySize);
    res = B_inv.array() * X.array();
    for(size_t i = 0; i < xSize; i++)
        for(size_t j = 0; j < ySize; j++){
            queueUnboundedMultRandom.push(rand(i, j));
            queueUnboundedMultRandom.push(res(i, j));
        }
    return;
}

void DoubleRandom::generate_mult_reduced_random_sharing_pairs(size_t num)
{
    DoubleShareBundle R(num, 2);
    R.reduced_random();
    gfpMatrix r = R.shares;
    
    ShareBundle prod(num, 1);
    prod.shares = r.col(0).array() * r.col(1).array();
    prod.reduce_degree();

    MultVerifier::push_triples(r.col(0), r.col(1), prod.shares, true);

    for(size_t i = 0; i < num; i++){
        queueMultReducedRandomPair.push(R.shares(i, 0));
        queueMultReducedRandomPair.push(R.aux_shares(i, 0));
        queueMultReducedRandomPair.push(R.shares(i, 1));
        queueMultReducedRandomPair.push(R.aux_shares(i, 1));
        queueMultReducedRandomPair.push(prod.shares(i, 0));
    }
    
    return;
}

void DoubleRandom::generate_reused_mult_reduced_random_sharing_pairs(size_t num)
{
    DoubleShareBundle R(num, 5);
    R.reduced_random();
    gfpMatrix r = R.shares;
    
    ShareBundle prod(num, 4);
    for (size_t i = 0; i < 4; i++)
    {
        prod.shares.col(i) = r.col(0).array() * r.col(i+1).array();
    }
    prod.reduce_degree();

    MultVerifier::push_reused_triples(r.col(0), r(seq(0, last), seqN(1, 4)), prod.shares, true);

    for(size_t i = 0; i < num; i++){
        gfpMatrix tmp(2, 7);
        tmp(0, seqN(0, 5)) = R.shares(i, seq(0, last));
        tmp(1, seqN(0, 5)) = R.aux_shares(i, seq(0, last));
        tmp(0, 5) = prod.shares(i, 0);
        tmp(1, 5) = prod.shares(i, 1);
        tmp(0, 6) = prod.shares(i, 2);
        tmp(1, 6) = prod.shares(i, 3);
        queueReusedMultReducedRandomPair.push(tmp);
    }
    
    return;
}

void DoubleRandom::generate_reused_mult_reduced_random_sharing_pairs_16(size_t num)
{
    DoubleShareBundle R(num, 17);
    R.reduced_random();
    gfpMatrix r = R.shares;
    
    ShareBundle prod(num, 16);
    for (size_t i = 0; i < 16; i++)
    {
        prod.shares.col(i) = r.col(0).array() * r.col(i+1).array();
    }
    prod.reduce_degree();

    MultVerifier::push_reused_triples(r.col(0), r(seq(0, last), seqN(1, 16)), prod.shares, true);

    for(size_t i = 0; i < num; i++){
        gfpMatrix tmp(2, 25);
        tmp(0, seqN(0, 17)) = R.shares(i, seq(0, last));
        tmp(1, seqN(0, 17)) = R.aux_shares(i, seq(0, last));
        for (size_t j = 0, k = 0; j < 16; j+=2, k++)
        {
            tmp(0, 17+k) = prod.shares(i, j);
            tmp(1, 17+k) = prod.shares(i, j+1);
        }
        queueReusedMultReducedRandomPair16.push(tmp);
    }
    
    return;
}

void DoubleRandom::generate_beaver_triple(size_t dim_x, size_t dim_y, size_t dim_z)
{
    // Generate the random sharings
    ShareBundle R(dim_x*dim_y+dim_y*dim_z, 1);
    R.random();

    // Compute the Matrix Multiplication
    gfpMatrix A = R.shares(seqN(0, dim_x*dim_y), 0).reshaped<Eigen::RowMajor>(dim_x, dim_y);
    gfpMatrix B = R.shares(seqN(dim_x*dim_y, dim_y*dim_z), 0).reshaped<Eigen::RowMajor>(dim_y, dim_z);
    ShareBundle CShareBundle(dim_x*dim_z, 1);
    CShareBundle.shares = (A * B).reshaped<Eigen::RowMajor>();
    CShareBundle.reduce_degree();
    MatrixTriple mat(A, B, CShareBundle.shares.reshaped<Eigen::RowMajor>(dim_x, dim_z));
    if(MultVerifier::verify_protocol.compare("Falcon") == 0) { // When use Falcon, it needs to provide correct triple
        MultVerifier::push_matrix_triple(mat.left, mat.right, mat.product, true);
    }
    queueBeaverTriples.push(mat);
}

void DoubleRandom::generate_beaver_triples(size_t num, vector<MatrixDimension> vectorDimensions)
{
    size_t random_num = 0;
    size_t C_num = 0;
    for(MatrixDimension dim : vectorDimensions) {
        random_num += dim.both_dimension();
        C_num += dim.z_dimension();
    }
    assert(random_num == num);

    // Generate the random sharings
    ShareBundle R(num, 1);
    R.random();

    // Compute the Matrix Multiplication
    size_t index = 0;
    size_t index_C = 0;
    ShareBundle CShareBundle(C_num, 1);
    vector<MatrixTriple> vectorMatrixTriples;
    for(MatrixDimension dim : vectorDimensions) {
        gfpMatrix A = R.shares(seqN(index, dim.dim_x*dim.dim_y), 0).reshaped<Eigen::RowMajor>(dim.dim_x, dim.dim_y);
        gfpMatrix B = R.shares(seqN(index+dim.dim_x*dim.dim_y, dim.dim_y*dim.dim_z), 0).reshaped<Eigen::RowMajor>(dim.dim_y, dim.dim_z);
        vectorMatrixTriples.push_back(MatrixTriple(A, B));
        index += dim.both_dimension();

        CShareBundle.shares(seqN(index_C, dim.dim_x*dim.dim_z), 0) = (A * B).reshaped<Eigen::RowMajor>();
        index_C += dim.z_dimension();
    }
    CShareBundle.reduce_degree();

    // Phase the product
    index_C = 0;
    for(int i = 0; i < vectorMatrixTriples.size(); i++) {
        vectorMatrixTriples[i].product = CShareBundle.shares(seqN(index_C, vectorDimensions[i].dim_x*vectorDimensions[i].dim_z), 0).reshaped<Eigen::RowMajor>(vectorDimensions[i].dim_x, vectorDimensions[i].dim_z);
        index_C += vectorDimensions[i].z_dimension();
        if(MultVerifier::verify_protocol.compare("Falcon") == 0) { // When use Falcon, it needs to provide correct triple
            MultVerifier::push_matrix_triple(vectorMatrixTriples[i].left, vectorMatrixTriples[i].right, vectorMatrixTriples[i].product, true);
        }
        queueBeaverTriples.push(vectorMatrixTriples[i]);
    }

    return;
}


}

