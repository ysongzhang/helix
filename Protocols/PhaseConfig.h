#ifndef PROTOCOLS_PHASE_CONFIG_H_
#define PROTOCOLS_PHASE_CONFIG_H_
#include "Networking/Player.h"
#include "Types/Triple.h"
namespace hmmpc
{

class PhaseConfig
{
    size_t offlineSent;
    mutable double offlineTimer;//
    mutable NamedCommStats offlineComm;

    size_t onlineSent;
    mutable double onlineTimer;//
    mutable NamedCommStats onlineComm;
    
    // * Note 
    // We'd generate the various random sharings in the true offline phase (totally before the online phase).
    // But it's trivial to manage the exact number of all the random sharings in some complicated operations, e.g. prefix_op
    // Hence, we'd like to delte the offline cost out of the online cost.
    bool trueOffline = false; 
    bool isOffline = false;
    bool isOnline = false;

    size_t cntRandom = 0;
    size_t cntRandomBit = 0;

    size_t cntReducedRandom = 0;
    size_t cntTruncatedRandom = 0;
    size_t cntTruncatedRandomInML = 0;
    size_t cntReducedTruncatedRandom = 0;
    size_t cntReducedTruncatedInMLRandom = 0;
    size_t cntUnboundedMultRandom = 0;
    size_t constUnboundedSize = 0;

    size_t cntRTRandomWithDifferentPrecision = 0;

    size_t cntMultReducedRandomPair = 0;
    size_t cntReusedMultReducedRandomPair = 0;
    size_t cntReusedMultReducedRandomPair16 = 0;

    size_t cntBeaverMultRandom = 0;
    queue<MatrixDimension> queueBeaverMultDimension;
    queue<MatrixDimension> queueBeaverMultDimensionOffline;
public:
    PhaseConfig():offlineSent(0), offlineTimer(0), onlineSent(0), onlineTimer(0){}
    void init(int n, int t, ThreadPlayer *_P, int _Pking = 0);
    void set_input_file(string fn);
    void close_input_file();

    bool is_true_offline(){return trueOffline;}
    void setTrueOffline(){trueOffline = true;}
    void offTrueOffline(){trueOffline = false;}
    bool is_Offline(){return isOffline;}
    bool is_Online(){return isOnline;}

    void start_offline();
    void end_offline();
    void print_offline_communication();
    void clear_offline_status();

    void switch_to_offline();
    void switch_to_online();

    void generate_random(string filename);// read argv from filename
    void generate_random_sharings(size_t n);
    void generate_random_bits(size_t n);
    void generate_reduced_random_sharings(size_t n);
    void generate_truncated_random_sharings(size_t n);
    void generate_truncated_random_sharings(size_t n, size_t precision);
    void generate_reduced_truncated_sharings(size_t n);
    void generate_reduced_truncated_sharings(size_t n, size_t logLearningRate, size_t logMiniBatch);// used to update parameters in ML
    void generate_reduced_truncated_sharings(size_t n, size_t precision);
    void generate_reduced_truncated_sharings(size_t n, vector<size_t> &precision, size_t num_repetition);

    void generate_unbounded_mult_random_sharings(size_t xSize, size_t ySize);
    void generate_unbounded_mult_random_sharings(size_t num);

    void generate_mult_reduced_random_sharing_pairs(size_t n);
    void generate_reused_mult_reduced_random_sharing_pairs(size_t n);
    void generate_reused_mult_reduced_random_sharing_pairs_16(size_t n);

    void generate_beaver_triple(size_t dim_x, size_t dim_y, size_t dim_z, bool is_offline = false);
    void generate_beaver_triples(size_t n, int triple_num, ifstream &in);

    void start_online();
    void end_online();
    void print_online_communication();
    void clear_online_status();

    size_t get_global_communication(string phase);

    void print_communication_oneline();
    void print_online_comm_online();
};

}

#endif