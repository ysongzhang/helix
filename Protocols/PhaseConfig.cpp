#include "PhaseConfig.h"
#include "Protocols/RandomShare.h"
namespace hmmpc
{

void PhaseConfig::init(int n, int t, ThreadPlayer *_P, int _Pking)
{
    ShareBase::n_players = n;
    ShareBase::threshold = t;
    ShareBase::P = _P;
    ShareBase::Pking = _Pking;
    ShareBase::Phase = this; // Bind the phase pointer.

    // Suppose this the the random seed agreed among the parties.
    const char key[SEED_SIZE]="ThisIs7heSeed.";
    ShareBase::PRNG_agreed.SetSeed((const octet*)key);
    // cout<<ShareBase::PRNG_agreed.get_uint()<<endl;

    ShareBase::init_vandermondes();
    ShareBase::init_reconstruction_vectors();

    ShareBase::init_bits_coeff();

    ShareBase::send_buffers.reset(n);
    ShareBase::receive_buffers.reset(n);

    ShareBase::send_buffers_PRG.reset(t+1);
    ShareBase::receive_buffers_PRG.reset(n);
}

void PhaseConfig::set_input_file(string fn)
{
    ShareBase::set_input_file(fn);
}

void PhaseConfig::close_input_file()
{
    ShareBase::close_input_file();
}

void PhaseConfig::start_offline()
{
    isOffline = true;
    auto &P = ShareBase::P;
    P->comm_stats.clear();
    P->timer.reset();
    P->sent = 0;
    P->timer.start();
}

void PhaseConfig::end_offline()
{
    isOffline = false;
    auto &P = ShareBase::P;
    P->timer.stop();
    offlineTimer += P->timer.elapsed();
    offlineSent += P->sent;
    offlineComm += P->comm_stats;
}

void PhaseConfig::clear_offline_status()
{
    offlineTimer = 0;
    offlineSent = 0;
    offlineComm.clear();
}

void PhaseConfig::print_offline_communication()
{
    auto &P = ShareBase::P;
    cout << "---OFFLINE----"<<endl;
    cout << "Time (for P"<<P->my_num()<<") = "<<offlineTimer<<" seconds"<<endl;
    cout << "Data sent (for P"<<P->my_num()<<") = "<< offlineSent/1e6 <<" MB in ~"<<offlineComm.total_rounds()<<" rounds"<<endl;
    octetStream o;
    o.store(offlineSent);
    P->send_all_no_stats(o);
    octetStreams os;
    P->receive_all_no_stats(os);
    size_t global = offlineSent;
    for(int i = 0; i < P->num_players(); i++)
        if(i != P->my_num())
            global += os[i].get_int(8);
    cout << "Global data sent = "<<global/1e6<<" MB (all parties)"<<endl;
    if(!is_true_offline()){
        cout <<"Preprocessed Randomness:"<<endl;

        // For offline files
        cout <<"#Random = "<<cntRandom<<endl;
        cout <<"#DoubleR = "<<cntReducedRandom<<endl;
        cout <<"#RandomBit = "<<cntRandomBit<<endl;
        cout <<"#RandomUn = "<<cntUnboundedMultRandom<<" "<<constUnboundedSize<<endl;
        cout <<"#RandomTrunc = "<<cntTruncatedRandom<<endl;
        cout <<"#RandomReTrunc = "<<cntReducedTruncatedRandom<<endl;
        cout <<"#MultReducedPair = "<<cntMultReducedRandomPair<<endl;
        cout <<"#ReMultReducedPair = "<<cntReusedMultReducedRandomPair<<endl;
        cout <<"#ReMultReducedPair16 = "<<cntReusedMultReducedRandomPair16<<endl;
        int triple_num = queueBeaverMultDimension.size() + queueBeaverMultDimensionOffline.size();
        cout <<"#BeaverTriplesNum = "<<triple_num<<" "<<cntBeaverMultRandom<<endl;
        int index = 0;
        for (int i = 0; i < queueBeaverMultDimensionOffline.size(); index++, i++)
        {
            MatrixDimension dim = queueBeaverMultDimensionOffline.front();
            cout <<"#Index "<<(i+1)<<" = "<<dim.dim_x<<" "<<dim.dim_y<<" "<<dim.dim_z<<" "<<endl;
            queueBeaverMultDimensionOffline.pop();
        }
        for (int i = index; i < triple_num; i++)
        {
            MatrixDimension dim = queueBeaverMultDimension.front();
            cout <<"#Index "<<(i+1)<<" = "<<dim.dim_x<<" "<<dim.dim_y<<" "<<dim.dim_z<<" "<<endl;
            queueBeaverMultDimension.pop();
        }
        
        // cout <<"#RandomTrunc(ML) = "<<cntTruncatedRandomInML<<endl;
        // cout <<"#RandomReTrunc(ML) = "<<cntReducedTruncatedInMLRandom<<endl;
        // cout <<"#RandomReTrunc(diff) = "<<cntRTRandomWithDifferentPrecision<<endl;
    }
    
    cout << "-------------"<<endl;
}

void PhaseConfig::start_online()
{
    isOnline = true;
    auto &P = ShareBase::P;
    P->comm_stats.clear();
    P->timer.reset();
    P->sent = 0;
    P->timer.start();
}



void PhaseConfig::end_online()
{
    isOnline = false;
    auto &P = ShareBase::P;
    P->timer.stop();
    onlineTimer += P->timer.elapsed();
    onlineComm += P->comm_stats;
    onlineSent += P->sent;
}

void PhaseConfig::clear_online_status()
{
    onlineTimer = 0;
    onlineComm.clear();
    onlineSent = 0;
}

void PhaseConfig::print_online_communication()
{
    auto &P = ShareBase::P;
    cout << "---ONLINE---"<<endl;
    cout << "Time (for P"<<P->my_num()<<") = "<<onlineTimer<<" seconds"<<endl;
    cout << "Data sent (for P"<<P->my_num()<<") = "<< onlineSent/1e6 <<" MB in ~"<<onlineComm.total_rounds()<<" rounds"<<endl;
    octetStream o;
    o.store(onlineSent);
    P->send_all_no_stats(o);
    octetStreams os;
    P->receive_all_no_stats(os);
    size_t global = onlineSent;
    for(int i = 0; i < P->num_players(); i++)
        if(i != P->my_num())
            global += os[i].get_int(8);
    cout << "Global data sent = "<<global/1e6<<" MB (all parties)"<<endl;
    cout << "-------------"<<endl;
}
size_t PhaseConfig::get_global_communication(string phase)
{
    size_t nSent;
    if(phase.compare("online")==0){
        nSent = onlineSent;
    }else if (phase.compare("offline")==0){
        nSent = offlineSent;
    }else if(phase.compare("all")==0){
        nSent = offlineSent+onlineSent;
    }else{
        assert(false && "get_global_communication");
    }

    auto &P = ShareBase::P;
    octetStream o;
    o.store(nSent);
    P->send_all_no_stats(o);
    octetStreams os;
    P->receive_all_no_stats(os);
    size_t global = nSent;
    for(int i = 0; i < P->num_players(); i++)
        if(i != P->my_num())
            global += os[i].get_int(8);
    return global;
}

void PhaseConfig::print_communication_oneline()
{
    auto &P = ShareBase::P;
    // communication per party
    cout<<onlineTimer<<" "<<get_global_communication("online")/1e6/P->num_players()<<" ";
    cout<<offlineTimer<<" "<<get_global_communication("offline")/1e6/P->num_players()<<" ";
    cout<<get_global_communication("all")/1e6/P->num_players()<<endl;
}


void PhaseConfig::print_online_comm_online()
{
    cout<<onlineTimer<<" "<<get_global_communication("online")/1e6<<endl;
}

void PhaseConfig::switch_to_offline()
{
    assert(isOnline);
    end_online();
    start_offline();
}

void PhaseConfig::switch_to_online()
{
    assert(isOffline);
    end_offline();
    start_online();
}

void PhaseConfig::generate_random(string filename)
{
    ifstream in(filename);
// #Random = 31861
// #DoubleR = 62087
// #RandomBit = 31842
// #RandomUn = 0 0
// #RandomTrunc = 0
// #RandomReTrunc = 266
// #MultReducedPair = 3840
// #ReMultReducedPair = 3584
// #ReMultReducedPair16 = 3584
// #BeaverTriplesNum = 0 0
    for(int i = 0; i < 10; i++){
        int x;
        int y;
        in>>x;
        if(i == 3 || i == 9) {
            in>>y;
        }
        if(x){
            switch (i)
            {
            case 0:
                generate_random_sharings(x);
                break;
            case 1:
                generate_reduced_random_sharings(x);
                break;
            case 2:
                generate_random_bits(x);
                break;
            case 3:
                generate_unbounded_mult_random_sharings(x, y);
                break;
            case 4:
                generate_truncated_random_sharings(x);
                break;
            case 5:
                generate_reduced_truncated_sharings(x);
                break;
            case 6:
                generate_mult_reduced_random_sharing_pairs(x);
                break;
            case 7:
                generate_reused_mult_reduced_random_sharing_pairs(x);
                break;
            case 8:
                generate_reused_mult_reduced_random_sharing_pairs_16(x);
                break;
            case 9:
                generate_beaver_triples(y, x, in);
                break;
            default:
                break;
            }
        }
        
    }
}

void PhaseConfig::generate_random_sharings(size_t n)
{
    cntRandom += n;
    // RandomShare::generate_random_sharings(n);
    RandomShare::generate_random_sharings_PRG_active(n);
}

void PhaseConfig::generate_random_bits(size_t n)
{
    cntRandomBit += n;
    RandomShare::generate_random_bits_active(n);
}

void PhaseConfig::generate_reduced_random_sharings(size_t n)
{
    cntReducedRandom += n;
    // DoubleRandom::generate_reduced_random_sharings(n);
    DoubleRandom::generate_reduced_random_sharings_PRG_active(n);
}

void PhaseConfig::generate_truncated_random_sharings(size_t n)
{
    cntTruncatedRandom += n;
    DoubleRandom::generate_truncated_random_sharings(n);
}

void PhaseConfig::generate_truncated_random_sharings(size_t n, size_t precision)
{
    cntTruncatedRandomInML += n;
    DoubleRandom::generate_truncated_random_sharings(n, precision);
}

void PhaseConfig::generate_reduced_truncated_sharings(size_t n)
{
    cntReducedTruncatedRandom += n;
    DoubleRandom::generate_reduced_truncated_random_sharings(n);
}

// Combine the reduced truncation with the opration *learning_rate/batch_size.
void PhaseConfig::generate_reduced_truncated_sharings(size_t n, size_t logLearningRate, size_t logMiniBatch)
{
    cntReducedTruncatedInMLRandom += n;
    DoubleRandom::generate_reduced_truncated_random_sharings(DoubleRandom::queueReducedTruncatedInML, n, FIXED_PRECISION+logLearningRate+logMiniBatch);
}

void PhaseConfig::generate_reduced_truncated_sharings(size_t n, size_t precision)
{
    cntRTRandomWithDifferentPrecision += n;
    DoubleRandom::generate_reduced_truncated_random_sharings(DoubleRandom::queueReducedTruncatedWithPrecisionRandom ,n, precision);
}

void PhaseConfig::generate_reduced_truncated_sharings(size_t n, vector<size_t> &precision, size_t num_repetition)
{
    cntRTRandomWithDifferentPrecision += n * num_repetition;
    DoubleRandom::generate_reduced_truncated_random_sharings(DoubleRandom::queueReducedTruncatedWithPrecisionRandom ,n, precision, num_repetition);
}

void PhaseConfig::generate_unbounded_mult_random_sharings(size_t xSize, size_t ySize)
{
    cntUnboundedMultRandom += xSize;
    if(constUnboundedSize==0){constUnboundedSize=ySize;}
    else if(ySize!=constUnboundedSize)assert(false && "Unbounded ySize is inconsistent");
    DoubleRandom::generate_unbounded_random_sharings_active(xSize, ySize);
}

void PhaseConfig::generate_unbounded_mult_random_sharings(size_t num)
{
    // cntUnboundedMultRandom += num;
    DoubleRandom::generate_unbounded_random_sharings_active(num);
}

void PhaseConfig::generate_mult_reduced_random_sharing_pairs(size_t n)
{
    cntMultReducedRandomPair += n;
    DoubleRandom::generate_mult_reduced_random_sharing_pairs(n);
}

void PhaseConfig::generate_reused_mult_reduced_random_sharing_pairs(size_t n)
{
    cntReusedMultReducedRandomPair += n;
    DoubleRandom::generate_reused_mult_reduced_random_sharing_pairs(n);
}

void PhaseConfig::generate_reused_mult_reduced_random_sharing_pairs_16(size_t n)
{
    cntReusedMultReducedRandomPair16 += n;
    DoubleRandom::generate_reused_mult_reduced_random_sharing_pairs_16(n);
}

void PhaseConfig::generate_beaver_triple(size_t dim_x, size_t dim_y, size_t dim_z, bool is_offline)
{
    MatrixDimension dim(dim_x, dim_y, dim_z);
    cntBeaverMultRandom += dim.both_dimension();
    if(is_offline) {
        queueBeaverMultDimensionOffline.push(dim);
    } else {
        queueBeaverMultDimension.push(dim);
    }
    DoubleRandom::generate_beaver_triple(dim_x, dim_y, dim_z);
}

void PhaseConfig::generate_beaver_triples(size_t n, int triple_num, ifstream &in)
{
    vector<MatrixDimension> dims;
    for(int i = 0; i < triple_num; i++) {
        int dim_x, dim_y, dim_z;
        in>>dim_x;
        in>>dim_y;
        in>>dim_z;
        dims.push_back(MatrixDimension(dim_x, dim_y, dim_z));
    }
    DoubleRandom::generate_beaver_triples(n, dims);
}

}