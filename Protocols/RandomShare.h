#ifndef PROTOCOLS_RANDOM_SHARE_H_
#define PROTOCOLS_RANDOM_SHARE_H_

#include "Types/Triple.h"
#include "Protocols/Share.h"
#include "Protocols/ShareBundle.h"

namespace hmmpc
{

class RandomShare: public ShareBase  // static class
{
protected:
    static size_t bundles(size_t num);

public:

    // Queues to store the preprocessed random sharings
    static queue<gfpScalar> queueRandom; // [r]_t
    static queue<gfpScalar> queueRandomBit;

    static void generate_random_sharings(size_t num); // Output into the queue
    static void generate_random_sharings_active(size_t num); // Output into the queue
    static void generate_random_bits(size_t num); //Output into the queue
    static void generate_random_bits_active(size_t num); // Output into the queue
    
    static void generate_random_sharings_PRG(size_t num); // Output into the queue
    static void generate_random_sharings_PRG_active(size_t num); // Output into the queue

    // Output directly, execute in the current phase
    static void generate_passive_random_sharings(gfpMatrix &res); // Execute in the current phase; semi-honest security
    static void generate_active_random_sharings(gfpMatrix &res);  // Execute in the current phase; malicious security
    static void generate_active_open_randoms(gfpMatrix &res);  // Execute in the current phase; malicious security
    static void generate_passive_random_sharings_PRG(gfpMatrix &res); // Execute in the current phase; semi-honest security
    static void generate_active_random_sharings_PRG(gfpMatrix &res);  // Execute in the current phase; malicious security
    static void generate_active_open_randoms_PRG(gfpMatrix &res);  // Execute in the current phase; malicious security
    // Output directly, execute in the online-offline paradigm
    static void generate_active_open_randoms_online_offline(gfpMatrix &res);  // Execute in the online-offline paradigm; malicious security

    
    static void get_random(queue<gfpScalar> &Q, gfpScalar &res);
    static void get_randoms(queue<gfpScalar> &Q, gfpMatrix &res);
    
};

inline size_t RandomShare::bundles (size_t num)
{
    size_t res = num / (n_players - threshold);
    if(num % (n_players-threshold)){res++;}
    return res;
}

class DoubleRandom: RandomShare  // static class
{
public:
    // Queue to store the preprocessed double sharings in order.
    static queue<gfpScalar> queueReducedRandom; // [r]_t, [r]_2t
    static queue<gfpScalar> queueTruncatedRandom; // [r/2^d]_t, [r]_t, [r_msb]_t
    static queue<gfpScalar> queueTruncatedRandomInML; // truncation: learning rate 2^-5 or 2^-7 ; batch size /2^7
    static queue<gfpScalar> queueReducedTruncatedRandom; // [r/2^d]_t, [r]_2t, [r]_t, [r_msb]_t
    static queue<gfpScalar> queueReducedTruncatedInML; // reduced + truncation: 2^-d; learning rate 2^-5 or 2^-7 ; batch size /2^7
    static queue<gfpScalar> queueUnboundedMultRandom; //([b1], [b1^-1]), ([bi], [bi-1 * bi^-1]) for i = 2, ..., l. 
    static queue<gfpScalar> queueMultReducedRandomPair; // ([r1]_t, [r1]_2t), ([r2]_t, [r2]_2t), [r1r2]_t
    static queue<gfpMatrix> queueReusedMultReducedRandomPair; // use matrix to store that contains 2 rows and 7 columns: ([r]_t, [r]_2t), ([r1]_t, [r1]_2t), ([r2]_t, [r2]_2t), ([r3]_t, [r3]_2t), ([r4]_t, [r4]_2t), [r*r1]_t, [r*r2]_t, [r*r3]_t, [r*r4]_t
    static queue<gfpMatrix> queueReusedMultReducedRandomPair16; // reused number is 16
    static queue<MatrixTriple> queueBeaverTriples; // [A]_t, [B]_t, [C]_t

    static queue<gfpScalar> queueReducedTruncatedWithPrecisionRandom; // For variable precision.

    // Output into the queue
    static void generate_reduced_random_sharings(size_t num); // [r]_t, [r]_2t
    static void generate_reduced_random_sharings_active(size_t num);
    static void generate_truncated_random_sharings(size_t num); // [r/2^d]_t, [r]_t, [r_msb]_t
    static void generate_truncated_random_sharings(size_t num, size_t precision); // [r/2^d]_t, [r]_t
    static void generate_reduced_truncated_random_sharings(size_t num); // [r/2^d]_t, [r]_2t, [r]_t, [r_msb]_t
    static void generate_reduced_truncated_random_sharings(queue<gfpScalar> &Q, size_t num, size_t precision); // [r/2^p]_t, [r]_2t, [r]_t
    static void generate_reduced_truncated_random_sharings(queue<gfpScalar> &Q, size_t num, vector<size_t> &precision, size_t num_repetitions);// different precision
    
    static void generate_unbounded_random_sharings(size_t num); // ([b1], [b1^-1]), ([bi], [bi-1 * bi^-1]) for i = 2, ..., l
    static void generate_unbounded_random_sharings_active(size_t num);
    static void generate_unbounded_random_sharings(size_t xSize, size_t ySize); // Each row corresponds an instance of unbounded prefix mult
    static void generate_unbounded_random_sharings_active(size_t xSize, size_t ySize);
    static void generate_non_zero_random_column_wise(gfpMatrix &res, gfpMatrix &prod_reveal, bool active = false);
    static void generate_non_zero_random_row_wise(gfpMatrix &res, gfpMatrix &prod_reveal, bool active = false);

    static void generate_reduced_random_sharings_PRG(size_t num);
    static void generate_reduced_random_sharings_PRG_active(size_t num);

    static void generate_mult_reduced_random_sharing_pairs(size_t num); // malicious security; ([r1]_t, [r1]_2t), ([r2]_t, [r2]_2t), [r1r2]_t
    static void generate_reused_mult_reduced_random_sharing_pairs(size_t num); // malicious security; ([r]_t, [r]_2t), ([r1]_t, [r1]_2t), ([r2]_t, [r2]_2t), ([r3]_t, [r3]_2t), ([r4]_t, [r4]_2t), [r*r1]_t, [r*r2]_t, [r*r3]_t, [r*r4]_t
    static void generate_reused_mult_reduced_random_sharing_pairs_16(size_t num); // malicious security;

    static void generate_beaver_triple(size_t dim_x, size_t dim_y, size_t dim_z); // malicious security;
    static void generate_beaver_triples(size_t num, vector<MatrixDimension> vectorDimensions);  // malicious security; The num is denoted as the number of random sharings which is needed.

    static void get_random_pair(queue<gfpScalar> &Q, gfpScalar &r, gfpScalar &aux_r);
    static void get_random_pairs(queue<gfpScalar> &Q, gfpMatrix &r, gfpMatrix &aux_r);
    
    static void get_random_triple(queue<gfpScalar>&Q, gfpScalar &r, gfpScalar &aux_r, gfpScalar &sub_r);
    static void get_random_triples(queue<gfpScalar> &Q, gfpMatrix &r, gfpMatrix &aux_r, gfpMatrix &sub_r);

    static void get_random_quadruple(queue<gfpScalar>&Q, gfpScalar &r, gfpScalar &aux_r, gfpScalar &aux_aux_r, gfpScalar &sub_r);
    static void get_random_quadruples(queue<gfpScalar> &Q, gfpMatrix &r, gfpMatrix &aux_r, gfpMatrix &aux_aux_r, gfpMatrix &sub_r);

    static void get_mult_reduced_random_sharing_pairs(gfpMatrix &r1, gfpMatrix &aux_r1, gfpMatrix &r2, gfpMatrix &aux_r2, gfpMatrix &prod);
    static void get_reused_mult_reduced_random_sharing_pairs(gfpMatrix &r1, gfpMatrix &aux_r1, gfpMatrix &r2, gfpMatrix &aux_r2, gfpMatrix &prod);
    static void get_reused_mult_reduced_random_sharing_pairs_16(gfpMatrix &r1, gfpMatrix &aux_r1, gfpMatrix &r2, gfpMatrix &aux_r2, gfpMatrix &prod);

    static void get_beaver_triple(gfpMatrix &a, gfpMatrix &b, gfpMatrix &c);
};

}
#endif