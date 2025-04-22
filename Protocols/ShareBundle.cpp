#include "Protocols/ShareBundle.h"
#include "Protocols/RandomShare.h"
#include "Protocols/Bit.h"
#include "Math/constMatrix.h"
#include "Protocols/BeaverTriper.h"
#include "Protocols/InputVerifier.h"
#include "Protocols/RevealVerifier.h"
#include "Protocols/MultVerifier.h"
#include "Types/maliciousCheck.h"
using Eigen::RowMajor;
using Eigen::seq;
using Eigen::seqN;
using Eigen::last;
namespace hmmpc
{

/*******************************************************************************************************
 * 
 *       Definition of static member and member functions about ShareBundle
 * The member functions of ShareBundle are almost similar to the Share that operate a bundle of sharings.
 * 
 * *****************************************************************************************************/

/*************************************************
 * 
 *       Basic operation about the sharings
 * 
 * ***********************************************/

/**
 * @brief Calculate and pack a bundle of sharings, which are of the same degree.
 * We firstly use the Vandermonde matrix to calcualte the sharings of the corresponding secrets.
 * Then we store my own shares into 'shares', which will be stored into the specific instance later.
 * As for other parties shares, we pack them into the octetStreams os[n], each of which corresponds to one party.
 * 
 * @param secrets [in]
 * @param degree [in]
 * @param shares [out]
 * @param os [out]
 */
void ShareBundle::calculate_sharings(const gfpMatrix &secrets, const int &degree, gfpMatrix &shares, octetStreams &os)
{
    gfpMatrix &vandermonde = (degree==threshold)? vandermonde_t : vandermonde_2t;
    gfpMatrix random_coeffs(degree, secrets.size());// Each column corresponds to the coefficients of one random polynomial of such degree. (Except the constant term)
    random_matrix(random_coeffs); 
    gfpMatrix sharings = vandermonde * random_coeffs;// Each column corresponds to one sharing of the secret
    for(size_t i = 0; i < secrets.size(); i++){
        sharings.col(i).array() += secrets(i);
    }
    // Extract my own shares and store them.
    shares = sharings.row( P->my_num()).reshaped<Eigen::AutoOrder>(shares.rows(), shares.cols());
    // TODO pack_row: add another functionality that not pack os[P->my_num()]
    pack_row(sharings, os);
}

/**
 * @brief Calculate the secrets from the corresponding sharings using the Lagrange Interpolation.
 * The reconstruction vector to point=0 has been evaluated in preprocessing stage.
 * 
 * @param secrets 
 * @param degree 
 * @param shares 
 * @param os 
 */
void ShareBundle::calculate_secrets(gfpMatrix &secrets, const int &degree, const gfpMatrix &shares, octetStreams &os)
{
    gfpVector &reconstruction = (degree==threshold)?reconstruction_vector_t:reconstruction_vector_2t;
    gfpMatrix sharings(degree+1, secrets.size());
    unpack_row(sharings, os);
    // If P is in the reconstruction set, the shares of P is not in os because we avoid the loopback.
    // Hence, we need to fill our own shares into the sharings.row(num)
    if(is_in_reconstruction_set(degree)){
        sharings.row(P->my_num()) = shares.reshaped<Eigen::RowMajor>().transpose();  // reshaped(): transfer to linear vector.
    }
    secrets = (reconstruction.transpose() * sharings).reshaped<Eigen::RowMajor>(secrets.rows(), secrets.cols());
}

void ShareBundle::distribute_sharings()
{
    octetStreams os(n_players);
    calculate_sharings(secrets, degree, shares, os);
    P->send_respective(os);
}

void ShareBundle::reconstruct_secrets()
{
    octetStreams os(degree+1);
    for(size_t player_no = 0; player_no < degree+1; player_no++)
    if(P->my_num()!=player_no){
        // P->receive_player(player_no, os[player_no]);
        P->request_receive(player_no, os[player_no]); // Requst to receive from player_no
    }

    for(size_t player_no = 0; player_no < degree+1; player_no++)
    if(P->my_num()!=player_no){
        P->wait_receive(player_no, os[player_no]);
    }

    calculate_secrets(secrets, degree, shares, os);
}

void ShareBundle::send_shares(int player_no)
{
    octetStream o;
    pack_matrix(shares, o);
    P->send_to(player_no, o);
}

void ShareBundle::receive_shares(int player_no)
{
    octetStream o;
    P->receive_player(player_no, o);
    unpack_matrix(shares, o);
}

void ShareBundle::send_shares_all()
{
    octetStream o;
    pack_matrix(shares, o);
    P->send_all(o);  // TODO_ZYS: 根据作者的DEBUG日志，似乎不能直接使用send_all函数，而是需要拆开request和wait来，后续如果通信出问题可以参考ShareBundle::reveal_dispersed()函数中的模式进行修改
}

void ShareBundle::receive_shares_all(octetStreams& os)
{
    P->receive_respective(os);
}

void ShareBundle::send_secrets()
{
    octetStream o;
    pack_matrix(secrets, o);
    P->send_all(o);
}

void ShareBundle::send_secrets(octetStream& o)
{
    pack_matrix(secrets, o);
    P->send_all(o);
}

void ShareBundle::receive_secrets(int player_no)
{
    octetStream o;
    P->receive_player(player_no, o);
    unpack_matrix(secrets, o);
}

void ShareBundle::receive_secrets(int player_no, octetStream& o)
{
    P->receive_player(player_no, o);
    unpack_matrix(secrets, o);
}

/**
 * @brief Calculate the block of secrets from the received octetStream.
 * The block starts from start_row, containing n_rows.
 * (block version of 'calculate_secrets')
 * 
 * @param degree 
 * @param start_row The sharings corresponds to the block of the matrix.
 * @param n_rows 
 * @param os The relevalent shares stored in the octetStreams. (size = degree + 1)
 */
void ShareBundle::calculate_block_secrets(const int &degree, const int &start_row, const int &n_rows, octetStreams &os)
{
    gfpVector &reconstruction = (degree==threshold)?reconstruction_vector_t:reconstruction_vector_2t;
    gfpMatrix sharings(degree+1, n_rows * cols());
    unpack_row(sharings, os);
    // Block starts from 'start_row', containing 'n_rows' rows.
    secrets(seqN(start_row, n_rows), seq(0, last)) = (reconstruction.transpose() * sharings).reshaped<RowMajor>(n_rows, cols()); 
}

/**
 * @brief Caculate the block of sharings, which starts from 'start_row' containing 'n_rows'.
 * (block version of 'calculate_sharings')
 * @param degree 
 * @param start_row 
 * @param n_rows 
 * @param os 
 */
void ShareBundle::calculate_block_sharings(const int &degree, const int &start_row, const int &n_rows, octetStreams &os)
{
    const gfpMatrix &vandermonde = (degree==threshold) ? vandermonde_t : vandermonde_2t;
    size_t n_secrets = n_rows * cols();
    gfpMatrix random_coeffs(degree, n_secrets);
    random_matrix(random_coeffs);
    gfpMatrix sharings = vandermonde * random_coeffs;
    
    // Secret ID in secrets matrix.
    size_t s_idx = start_row * cols();
    for(size_t i = 0; i < n_secrets; i++, s_idx++){
        sharings.col(i).array() += secrets(s_idx);
    }
    // Extract my own shares and store them into the block, which start from start_row containing n_rows.
    shares(seqN(start_row, n_rows), seq(0, last)) = sharings.row(P->my_num()).reshaped<RowMajor>(n_rows, cols());
    pack_row(sharings, os);
}

// * With PRG

/**
 * @brief Pi distribute a 2t-sharing. (ShareBundle Version of Share::distribute_2t_sharing_PRG)
 * Pi generates 2t shares from PRG. Pi calculate its own share from the 2t shares and its secret.
 * The communication is 0.
 *
 * Other parties can get its share directly from the PRG.
 * @param _secrets 
 * @param _shares 
 */
void ShareBundle::calculate_2t_sharings_PRG(const gfpMatrix&_secrets, gfpMatrix &_shares)
{
    int _degree = threshold<<1;
    gfpMatrix shares_prng(_degree, _secrets.size());
    random_matrix(shares_prng, PRNG_agreed);

    gfpMatrix material(_degree+1, _secrets.size());
    material.row(0) = _secrets.reshaped<Eigen::RowMajor>().transpose();
    material.bottomRows(_degree) = shares_prng;

    _shares = ( reconstruction_with_secret_2t.row(P->my_num()) * material).
                reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
}

/**
 * @brief Pi distribute a t-sharing. (ShareBundle Version of Share::distribute_t_sharing_PRG)
 * Pi generates t shares from PRG. We stipulcate that the t shares are P1 ... Pt.
 * Pi calculates t+1 shares for Pt+1 ... P2t, and P0, from the t shares from PRG and its secret.
 * 
 * P1 ... Pt can get its share from PRG without communication.
 * Pt+1 ... P2t, and P0 need receive from Pi. (Assume Pi is not in them)
 * 
 * @param _secrets 
 * @param _shares 
 * @param os stores the shares of Pt+1, ..., P2t, and P0.
 */
void ShareBundle::calculate_t_sharings_PRG(const gfpMatrix&_secrets, gfpMatrix &_shares, octetStreams &os)
{
    assert(degree==threshold);

    // Mapped party
    // P0 P1 ... Pt Pt+1

    // shares from prng: P2 ... Pt+1
    gfpMatrix shares_prng(degree, _secrets.size());
    random_matrix(shares_prng, PRNG_agreed);

    // material: P0 P2 ... Pt+1
    // True idx is P1 ... Pt
    gfpMatrix material(n_party_PRG(), _secrets.size());
    material.row(0) = _secrets.reshaped<Eigen::RowMajor>().transpose();
    material.bottomRows(degree) = shares_prng;

    // sharing(t+1, size()): Calculate shares for P_t+2 ... P_2t+1 P1
    // True idx is P_{t+1} ... P_{2t} P0
    gfpMatrix sharings = reconstruction_with_secret_t * material;

    if(P->my_num()>=1 && P->my_num()<=degree){
        // Obtain shares directly from PRG
        _shares = shares_prng.row(P->my_num()-1).reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
    }else if(P->my_num()==0){
        // P0 can get from the sharings that calculated.
        _shares = sharings.row(degree).reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
    }else{
        _shares = sharings.row(P->my_num() - degree - 1).reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
    }

    // Pack the sharings that contains shares of P_{t+1} ... P_{2t} P0.
    assert(sharings.rows()==os.size());
    pack_rows(sharings, os);
}

void ShareBundle::get_2t_sharings_PRG(int player_no, gfpMatrix &_shares)
{
    gfpMatrix shares_prng(threshold<<1, _shares.size());
    random_matrix(shares_prng, PRNG_agreed);

    if(P->my_num()>player_no){
        _shares = shares_prng.row(P->my_num() - 1).reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
    }else{
        _shares = shares_prng.row(P->my_num()).reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
    }
}

void ShareBundle::get_t_sharings_PRG(int player_no, gfpMatrix &_shares)
{
    gfpMatrix shares_prng(degree, _shares.size());
    random_matrix(shares_prng, PRNG_agreed);

    if(P->my_num()>=1 && P->my_num()<=degree){
        _shares = shares_prng.row(P->my_num()-1).reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
    }else{
        octetStream o;
        P->receive_player(player_no, o);
        unpack_matrix(_shares, o);
    }
    return;
}

void ShareBundle::distribute_sharings_PRG(const gfpMatrix &_secrets, gfpMatrix &_shares)
{
    if(degree==threshold){
        octetStreams os(degree+1);
        calculate_t_sharings_PRG(_secrets, _shares, os);
        P->send_respective(start_party_PRG(), n_party_PRG(), os);
    }else{
        calculate_2t_sharings_PRG(_secrets, _shares);
    }
}

void ShareBundle::get_sharings_PRG(int player_no, gfpMatrix &_shares)
{
    if(degree==threshold){get_t_sharings_PRG(player_no, _shares);}
    else{get_2t_sharings_PRG(player_no, _shares);}
}

// We divide get_t_sharings_PRG into two phase: the request phase and the wait phase.
void ShareBundle::get_t_sharings_PRG_request(int player_no, gfpMatrix &shares_prng, octetStream &o_receive)
{
    assert(shares_prng.rows()==degree);
    assert(shares_prng.cols()==size());
    random_matrix(shares_prng, PRNG_agreed);

    if(P->my_num()>=1 && P->my_num()<=degree){
        return;
    }else{
        P->request_receive(player_no, o_receive);
    }
    return;
}

void ShareBundle::get_t_sharings_PRG_wait(int player_no, const gfpMatrix &shares_prng, gfpMatrix &_shares, octetStream &o_receive)
{
    if(P->my_num()>=1 && P->my_num()<=degree){
        _shares = shares_prng.row(P->my_num()-1).reshaped<Eigen::RowMajor>(_shares.rows(), _shares.cols());
        return;
    }else{
        P->wait_receive(player_no, o_receive);
        unpack_matrix(_shares, o_receive);
    }
    return;
}

// Row-grained of calculate_t_sharings_PRG
void ShareBundle::calculate_block_t_sharings_PRG(const int &start_row, const int &n_rows, octetStreams &os_send)
{
    assert(degree==threshold);

    size_t n_secrets = n_rows * cols();
    gfpMatrix shares_prng(degree, n_secrets);
    random_matrix(shares_prng, PRNG_agreed);

    gfpMatrix material(n_party_PRG(), n_secrets);
    material.row(0) = secrets.middleRows(start_row, n_rows).reshaped<Eigen::RowMajor>().transpose();
    material.bottomRows(degree) = shares_prng;

    gfpMatrix sharings = reconstruction_with_secret_t * material;  // size: n-t * n_secrets

    if(P->my_num()>=1 && P->my_num()<=degree){  // f(2) ... f(t+1)
        // Obtain shares directly from PRG
        shares.middleRows(start_row, n_rows) = shares_prng.row(P->my_num()-1).reshaped<Eigen::RowMajor>(n_rows, cols());
    }else if(P->my_num()==0){  // f(1)
        // P0 can get from the sharings that calculated.
        shares.middleRows(start_row, n_rows) = sharings.row(degree).reshaped<Eigen::RowMajor>(n_rows, cols());
    }else{
        shares.middleRows(start_row, n_rows) = sharings.row(P->my_num() - degree - 1).reshaped<Eigen::RowMajor>(n_rows, cols());
    }

    // Pack the sharings that contains shares of P_{t+1} ... P_{2t} P0.
    assert(sharings.rows()==os_send.size());
    pack_rows(sharings, os_send);
}

// Row-grained of get_t_sharings_PRG_request
void ShareBundle::get_block_t_sharings_PRG_request(int player_no, const int &start_row, const int &n_rows, gfpMatrix &shares_prng, octetStream &o_receive)
{
    size_t n_secrets = n_rows * cols();

    assert(shares_prng.rows()==degree);
    assert(shares_prng.cols()==n_secrets);
    random_matrix(shares_prng, PRNG_agreed);

    if(P->my_num()>=1 && P->my_num()<=degree){
        return;
    }else{
        P->request_receive(player_no, o_receive);
    }
    return;
}

// Row-grained of get_t_sharings_PRG_wait
void ShareBundle::get_block_t_sharings_PRG_wait(int player_no, const int &start_row, const int &n_rows, const gfpMatrix &shares_prng, octetStream &o_receive)
{
    if(P->my_num()>=1 && P->my_num()<=degree){
        shares.middleRows(start_row, n_rows) = shares_prng.row(P->my_num()-1).reshaped<Eigen::RowMajor>(n_rows, cols());
        return;
    }else{
        P->wait_receive(player_no, o_receive);
        unpack_rows(shares, start_row, n_rows, o_receive);
    }
    return;
}

/*************************************************************************************
 *       Main operation about the sharings
 * (We use the condition 'DISPERSE_PKING' to change the normal version to dispersed version.)
 * 1. reduce_degree(): reduce the sharings of degree 2t to degree t.
 * 2. truncate(): truncate the last d bits of the sharings of degree t.
 * 3. reduce_truncate(): truncate the last d bits of the sharings of degree 2t to degree t.
 * ***********************************************************************************/

/**
 * @brief The reduction protocol / multiplication protocol of the ShareBundle version
 * can reduce a bundle of sharings of degree 2t to a bundle of sharing of degree t.
 * 
 * @return ShareBundle& 
 */
ShareBundle& ShareBundle::reduce_degree()
{
    // Normal
    double_degree();
    assert(degree == threshold<<1);

    DoubleShareBundle R(rows(), cols());
    R.reduced_random(); // Get a bundle of reduced random sharings.
    shares += R.aux_shares;

    // * Original DN Protocol
    reveal(true);
    degree >>= 1;
    shares = secrets - R.shares;

    // * Pking collects the shares of x+r and reveal its t-sharing.
    // reveal_blocks_dispersed();
    // degree >>= 1;
    // input_blocks_dispersed_PRG();  // obtain the t-sharing shares

    // TODO_ZYS: Input correctness check?

    // shares = shares - R.shares;
    return *this;
}

/**
 * @brief The truncation protocol protocol of the ShareBundle version
 * can truncate the last d bits of the sharings of degree t.
 * 
 * @return ShareBundle& 
 */
ShareBundle& ShareBundle::truncate()
{
    assert(degree == threshold);
    DoubleShareBundle R(rows(), cols());
    ShareBundle r_msb(rows(), cols());
    R.truncated_random(r_msb);

    // Encode to make sure MSB(a)=0
    shares.array() += R.aux_shares.array() + ConstEncode;
    
    reveal_truncate(FIXED_PRECISION, true);
    
    ShareBundle is_overflow(rows(), cols());
    getMSB_matrix(secrets, is_overflow.shares);
    is_overflow.shares = ((gfpScalar)1 - r_msb.shares.array()) * is_overflow.shares.array();

    // Fix the truncation result and decode it.
    shares = secrets.array() - R.shares.array() + ConstGapInTruncation * is_overflow.shares.array() - ConstDecode;
    return *this;
}

// TODO: Truncate any precision with MSB(a)=0 (Need to modify)
ShareBundle& ShareBundle::truncate(size_t precision)
{
    assert(degree == threshold);
    DoubleShareBundle R(rows(), cols());
    R.truncated_random(precision);
    shares += R.aux_shares;
    
    reveal_truncate(precision, true);
    shares = secrets - R.shares;
    return *this;
}

/**
 * @brief This protocol of the ShareBundle version
 * can truncate the last d bits of the sharings of degree 2t to the sharings of degree t.
 * 
 * @return ShareBundle& 
 */
ShareBundle& ShareBundle::reduce_truncate()
{
    // Normal
    double_degree();
    assert(degree == threshold<<1);
    DoubleShareBundle R(rows(), cols());
    ShareBundle r_msb(rows(), cols());
    R.reduced_truncated_random(r_msb);

    shares.array() += R.aux_shares.array() + ConstEncode;

    reveal(true); // semi-honest; 2t-sharing
    degree>>=1;
    mult_checked.array() = secrets.array() - R.aux_aux_shares.array() - ConstEncode;  // t-sharing
    truncate_matrix(secrets, FIXED_PRECISION);

    ShareBundle is_overflow(rows(), cols());
    getMSB_matrix(secrets, is_overflow.shares);
    is_overflow.shares = (1 - r_msb.shares.array()) * is_overflow.shares.array();

    // Fix the truncation result.
    shares = secrets.array() - R.shares.array() + ConstGapInTruncation * is_overflow.shares.array() - ConstDecode;
    return *this;
}

// TODO: Need to modify it
ShareBundle& ShareBundle::reduce_truncate(size_t precision)
{
    double_degree();
    assert(degree == threshold<<1);
    DoubleShareBundle R(rows(), cols());
    R.reduced_truncated_random(precision);

    shares += R.aux_shares;
    reveal(true); // semi-honest; 2t-sharing
    degree>>=1;
    mult_checked = secrets - R.aux_aux_shares; // t-sharing
    truncate_matrix(secrets, precision);

    shares = secrets - R.shares;
    return *this;
}

// TODO: Need to modify it
// Note precision is as the same size with sharebundles.
ShareBundle& ShareBundle::reduce_truncate(vector<size_t> &precision)
{
    double_degree();
    assert(degree == threshold<<1);
    DoubleShareBundle R(rows(), cols());
    R.reduced_truncated_random(precision);

    shares += R.aux_shares;
    reveal(true); // semi-honest; 2t-sharing
    degree>>=1;
    mult_checked = secrets - R.aux_aux_shares; // t-sharing
    for(size_t i = 0; i < size(); i++){
        secrets(i).truncate(precision[i]);
    }

    shares = secrets - R.shares;
    return *this;
}

// Test: Two layer multiplication to compute [x]_2t * [y]_t
// First layer: compute ?*? = [x]_2t -> reduce degree [x]_t
// Second layer: compute [x]_t * [y]_t
ShareBundle& ShareBundle::reduce_degree_1stLayer(const ShareBundle &y, BeaverTriple &triple)
{
    ShareBundle* x = this;  // 2t-sharing

    assert(rows()==y.rows());
    assert(cols()==y.cols());
    assert(rows()==triple.rows());
    assert(cols()==triple.cols());

    // [rX] for [x]_2t and [rY] for [r]_t*[y]_t  
    DoubleShareBundle R(rows()<<1, cols());
    R.reduced_random();

    triple.a_share() = R.shares.topRows(rows()); // [a]_t = [rX]_t s.t. [x] = u - [rX]_t 
    triple.b_share() = - y.shares;               // [b]_t = [-y]_t s.t. [y] = 0 - (-[y]_t)  and v = 0
    triple.v_value().setConstant(0);

    ShareBundle combine(rows()<<1, cols());
    combine.shares.topRows(rows()) = shares + R.aux_shares.topRows(rows()); // [x]_2t + [rX]_2t, [u]_2t in paper
    combine.shares.bottomRows(rows()) = triple.a_share().array() * triple.b_share().array() // [rX]_t * [-y]_t = [- rX * y]_2t
                                        + R.aux_shares.bottomRows(rows()).array(); // + random of degree-2t; [u_i]_2t in paper
    
    
    combine.double_degree();//BUG LOG: The degree is 2t.
    combine.reveal(true);
    shares = combine.secrets.topRows(rows()) - R.shares.topRows(rows()); // [x] = e - [rX]_t
    triple.u_value() = combine.secrets.topRows(rows()); // u = e

    triple.c_share() = combine.secrets.bottomRows(rows()) - R.shares.bottomRows(rows()); // [c] = [ab]

    // Verify
    MultVerifier::push_triples(triple.a_share(), triple.b_share(), triple.c_share());

    return *this;
}

/**
 * @brief Generalize two layer-mult: compute ([xy1], [xy2], ...) from input [x]_2t and ([y1], [y2], ...)
 * [x]_2t is one input and [y_i] is the other input.
 * 
 * This method calculates two things.
 * 1. Reduce degree to finish the first-layer mult gate.
 *      obtain [x] = u - [rX] where [a] = [rX]
 * 2. Calculate the Beaver Triples for the second-layer mult gate
 *      for each y_i
 *      [y_i]_t = 0 - (-[y_i]) where v = 0 and [b] = -[y_i]
 *      compute the Beaver Triple ([a], [b], [c]) where [c] = [ab] = [rX * (-y_i)].
 *      So the evaluation of [c] can be done along with step 1.
 * @param y :
 * @param triples 
 * @return ShareBundle& 
 */
ShareBundle& ShareBundle::reduce_degree_1stLayer(const vector<ShareBundle>&y, vector<BeaverTriple>&triples)
{
    // ShareBundle* x = this;
    size_t len = y.size();
    assert(triples.size()==len);

    // The first 'rows()' double shares for [x]_2t, and the remaining for the evaluation of BeaverTriples.
    // [rX] = R.shares.topRows(rows());
    DoubleShareBundle R(rows()+rows()*len, cols());
    R.reduced_random();

    // Combine [x]_2t and [rX]*[-y_i]
    ShareBundle combine(rows()+rows()*len, cols());

    // [x+r]_2t = [x]_2t + [r]_2t
    combine.shares.topRows(rows()) = shares + R.aux_shares.topRows(rows());

    for(size_t i = 0, idxRow = rows(); i < len; i++, idxRow+=rows()){
        assert(rows()==y[i].rows());
        assert(cols()==y[i].cols());
        assert(rows()==triples[i].rows());
        assert(cols()==triples[i].cols());
        triples[i].a_share() = R.shares.topRows(rows());    // [x] = u - [rX] where [a] = [rX]
        triples[i].b_share() = -y[i].shares; // [y_i]_t = 0 - (-[y_i]) where v = 0 and [b] = -[y_i]
        triples[i].v_value().setConstant(0);

        combine.shares.middleRows(idxRow, rows()) = triples[i].a_share().array() * triples[i].b_share().array()
                                                    + R.aux_shares.middleRows(idxRow, rows()).array();
    }

    combine.double_degree();//BUG LOG: The degree is 2t.

    combine.reveal(true);
    
    // Step 1: obtain [x]_t = e - [rX]_t
    shares = combine.secrets.topRows(rows()) - R.shares.topRows(rows());
    // Step 2: obtain each BeaverTriple [c]_t = e_i - [rY_i]_t
    gfpMatrix checkedA(rows()*len, cols());
    gfpMatrix checkedB(rows()*len, cols());
    gfpMatrix checkedC(rows()*len, cols());
    for(size_t i = 0, idxRow = rows(); i < len; i++, idxRow+=rows()){
        triples[i].u_value() = combine.secrets.topRows(rows());
        triples[i].c_share() = combine.secrets.middleRows(idxRow, rows()) - R.shares.middleRows(idxRow, rows());
        checkedA(seqN(idxRow-rows(), rows()), seq(0, last)) = triples[i].a_share();
        checkedB(seqN(idxRow-rows(), rows()), seq(0, last)) = triples[i].b_share();
        checkedC(seqN(idxRow-rows(), rows()), seq(0, last)) = triples[i].c_share();
    }

    // Verify
    MultVerifier::push_triples(checkedA, checkedB, checkedC);
    return *this;
}


// TODO: Need to modify it
// It is used to update the parameters in ML.
// w = w - lr*X*delta/batch_size -> w - X*delta* (lr/batch_size)
ShareBundle& ShareBundle::reduce_truncate(size_t logLearningRate, size_t logMiniBatch)
{
    double_degree();
    assert(degree == threshold<<1);
    DoubleShareBundle R(rows(), cols());
    R.reduced_truncated_random(logLearningRate, logMiniBatch);

    shares += R.aux_shares;
    reveal(true); // semi-honest; 2t-sharing
    degree>>=1;
    mult_checked = secrets - R.aux_aux_shares; // t-sharing
    truncate_matrix(secrets, FIXED_PRECISION+logLearningRate+logMiniBatch);

    shares = secrets - R.shares;
    return *this;
}

/*********************************************************************
 * 
 *       Input Methods
 * 
 * There are serveral input methods with tiny difference.
 * The only difference is the source of the secret.
 * 
 * 1. input_from_file: secrets from the file
 * 2. input_from_party: secrets from the assigned value
 * 3. input_from_pking: the specific case of input_from_party
 * 4. input_from_random: secrets from the PRNG
 * 
 * *******************************************************************/

void ShareBundle::input_from_file(int player_no, bool active)
{
    if(player_no == P->my_num()){
        for(size_t i = 0; i < secrets.size(); i++){
            in>>secrets(i);
        }
        distribute_sharings();
    }else{receive_shares(player_no);}

    if (active && (degree==threshold)) {  // malicious security
        // Correctness check
        InputVerifier::push_matrix(shares);
    }   
}

void ShareBundle::input_from_party(int player_no, bool active)
{
    if(player_no == P->my_num()){distribute_sharings();}
    else { receive_shares(player_no);}

    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}

void ShareBundle::input_from_pking(bool active)
{
    input_from_party(Pking, active);
}

void ShareBundle::input_from_random(int player_no, bool active)
{
    if(player_no == P->my_num()){
        random_matrix(secrets);
        distribute_sharings();
    }else {receive_shares(player_no);}

    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}

/*********************************************************************
 * 
 *       Input Methods of multi-thread version
 * 
 * There are corresponing input methods of the multi-thread version.
 * The distributer calculates the sharings and only requests to send them.
 * The receiver only requests to receive the shares.
 * 
 * 1. input_from_file_request: request to input secrets from the file
 * 2. input_from_party_request: request to input secrets from the assigned value
 * 3. input_from_pking_request: the specific case of input_from_party
 * 4. input_from_random_request: request to input secrets from the PRNG
 * 
 * *******************************************************************/
void ShareBundle::input_from_file_request(int player_no, octetStreams &os_send, octetStream &o_receive)
{
    if(player_no == P->my_num()){
        for(size_t i = 0; i < secrets.size(); i++){
            in>>secrets(i);
        }
        calculate_sharings(secrets, degree, shares, os_send);
        P->request_send_respective(os_send);
    }else{P->request_receive(player_no, o_receive);}
}

// We use the send_buffers to distribute the sharings.
// We use to receive_buffers to receive the shares.
void ShareBundle::input_from_file_request(int player_no)
{
    // Clear the buffer.
    if(player_no == P->my_num()){send_buffers.clear();}
    else{receive_buffers[player_no].clear();}

    input_from_file_request(player_no, send_buffers, receive_buffers[player_no]);
}

/**
 * @brief The request part the input_from_random in the multi-thread version.
 * The distributer calculates the sharings and only requests to send them.
 * The receiver only requests to receive the shares.
 * @param player_no 
 * @param os_send 
 * @param o_receive 
 */
void ShareBundle::input_from_random_request(int player_no, octetStreams &os_send, octetStream &o_receive)
{
    if(player_no == P->my_num()){
        random_matrix(secrets);
        calculate_sharings(secrets, degree, shares, os_send);
        P->request_send_respective(os_send);
    }else{P->request_receive(player_no, o_receive);}
}

// We use the send_buffers to distribute the sharings.
// We use to receive_buffers to receive the shares.
void ShareBundle::input_from_random_request(int player_no)
{
    if(player_no == P->my_num()){send_buffers.clear();}
    else{receive_buffers[player_no].clear();}
    input_from_random_request(player_no, send_buffers, receive_buffers[player_no]);
}

void ShareBundle::input_from_party_request(int player_no, octetStreams &os_send, octetStream &o_receive)
{
    if(player_no == P->my_num()){
        calculate_sharings(secrets, degree, shares, os_send);
        P->request_send_respective(os_send);
    }else{P->request_receive(player_no, o_receive);}
}

void ShareBundle::input_from_party_request(int player_no)
{
    // Clear the buffer.
    if(player_no == P->my_num()){send_buffers.clear();}
    else{receive_buffers[player_no].clear();}

    input_from_party_request(player_no, send_buffers, receive_buffers[player_no]);
}

// When inputting from Pking,
// We use the send_buffers to distribute the sharings.
// We use to receive_buffers to receive the shares.
void ShareBundle::input_from_pking_request()
{
    // Clear the buffer.
    if(Pking == P->my_num()){send_buffers.clear();}
    else{receive_buffers[Pking].clear();}

    input_from_party_request(Pking, send_buffers, receive_buffers[Pking]);
}

/**
 * @brief The wait part of the multi-thread version for inputting sharings.
 * The distribtuer waits to send the shares to the corresponding party.
 * The receiver waits to receive the shares from the distributer.
 * @param player_no 
 * @param os_send 
 * @param o_receive 
 */
void ShareBundle::finish_input_from(int player_no, octetStreams &os_send, octetStream &o_receive, bool active)
{
    if(player_no == P->my_num()){
        P->wait_send_respective(os_send);
    }else{
        P->wait_receive(player_no, o_receive);
        unpack_matrix(shares, o_receive);
    }

    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}

// The wait part when inputting sharings,
// We use the send_buffers to distribute the sharings.
// We use to receive_buffers to receive the shares.
void ShareBundle::finish_input_from(int player_no, bool active)
{
    finish_input_from(player_no, send_buffers, receive_buffers[player_no], active);
}

// The wait part when inputting from Pking,
// We use the send_buffers to distribute the sharings.
// We use to receive_buffers to receive the shares.
void ShareBundle::finish_input_from_pking(bool active)
{
    finish_input_from(Pking, send_buffers, receive_buffers[Pking], active);
}

// With PRG
void ShareBundle::input_from_party_PRG(int player_no, bool active)
{
    if(player_no == P->my_num()){distribute_sharings_PRG(secrets, shares);}
    else{get_sharings_PRG(player_no, shares);}
    
    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}

void ShareBundle::input_from_random_PRG(int player_no, bool active)
{
    if(player_no == P->my_num()){
        random_matrix(secrets);
        distribute_sharings_PRG(secrets, shares);
    }else{
        get_sharings_PRG(player_no, shares);
    }

    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}

void ShareBundle::input_from_party_request_PRG(int player_no, octetStreams &os_send, gfpMatrix &shares_prng, octetStream &o_receive)
{
    assert(degree==threshold);
    if(player_no == P->my_num()){
        os_send.reset(n_party_PRG());
        // shares_prng.resize(degree, size());

        calculate_t_sharings_PRG(secrets, shares, os_send);
        P->request_send_respective(start_party_PRG(), n_party_PRG(), os_send);
    }else{
        shares_prng.resize(degree, size());
        get_t_sharings_PRG_request(player_no, shares_prng, o_receive);
    }
}

void ShareBundle::input_from_random_request_PRG(int player_no, octetStreams &os_send, gfpMatrix &shares_prng, octetStream &o_receive)
{
    assert(degree==threshold);
    if(player_no==P->my_num()){
        random_matrix(secrets);
        os_send.reset(n_party_PRG());
        // shares_prng.resize(degree, size());

        calculate_t_sharings_PRG(secrets, shares, os_send);
        P->request_send_respective(start_party_PRG(), n_party_PRG(), os_send);
    }else{
        shares_prng.resize(degree, size());
        get_t_sharings_PRG_request(player_no, shares_prng, o_receive);
    }
}

void ShareBundle::finish_input_from_PRG(int player_no, octetStreams &os_send, const gfpMatrix &shares_prng, octetStream &o_receive, bool active)
{
    assert(degree==threshold);
    if(player_no==P->my_num()){
        P->wait_send_respective(start_party_PRG(), n_party_PRG(), os_send);
    }
    else{
        get_t_sharings_PRG_wait(player_no, shares_prng, shares, o_receive);
    }

    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}

void ShareBundle::input_from_party_request_PRG(int player_no)
{
    assert(degree==threshold);
    if(player_no==P->my_num()){send_buffers_PRG.clear();}
    else{
        shares_buffers_PRG.resize(degree, size());
        receive_buffers_PRG[player_no].clear();
    }

    input_from_party_request_PRG(player_no, send_buffers_PRG, shares_buffers_PRG, receive_buffers_PRG[player_no]);
}

void ShareBundle::input_from_random_request_PRG(int player_no)
{
    assert(degree==threshold);
    if(player_no==P->my_num()){send_buffers_PRG.clear();}
    else{
        shares_buffers_PRG.resize(degree, size());
        receive_buffers_PRG[player_no].clear();
    }

    input_from_random_request_PRG(player_no, send_buffers_PRG, shares_buffers_PRG, receive_buffers_PRG[player_no]);
}

void ShareBundle::finish_input_from_PRG(int player_no, bool active)
{
    finish_input_from_PRG(player_no, send_buffers_PRG, shares_buffers_PRG, receive_buffers_PRG[player_no], active);
}

// Row-grained of input_from_party_request_PRG
void ShareBundle::input_block_from_party_request_PRG(int player_no, const int &startRow, const int &nRows, octetStreams &os_send, gfpMatrix &shares_prng, octetStream &o_receive)
{
    assert(degree==threshold);
    if(player_no==P->my_num()){
        os_send.reset(n_party_PRG());

        calculate_block_t_sharings_PRG(startRow, nRows, os_send);  // os_send stores the shares of Pt+1, ..., P2t, and P0.
        P->request_send_respective(start_party_PRG(), n_party_PRG(), os_send);  // send respectively according to the parameter start_party_PRG().
    }else{
        shares_prng.resize(degree, nRows*cols());
        get_block_t_sharings_PRG_request(player_no, startRow, nRows, shares_prng, o_receive);
    }
}

void ShareBundle::finish_input_block_from_party_PRG(int player_no, const int &startRow, const int &nRows, octetStreams &os_send, const gfpMatrix &shares_prng, octetStream &o_receive)
{
    assert(degree==threshold);
    if(player_no==P->my_num()){
        P->wait_send_respective(start_party_PRG(), n_party_PRG(), os_send);
    }
    else{
        get_block_t_sharings_PRG_wait(player_no, startRow, nRows, shares_prng, o_receive);
    }

    // Note: This finish method is only used for some rows of input (block).
}

// PRG Version of input_blocks_dispersed.
void ShareBundle::input_blocks_dispersed_PRG(bool active)
{
    size_t blk_rows = rows() / n_players;
    size_t first_blk_rows = rows() - blk_rows * (n_players - 1);

    octetStreams os_send(n_party_PRG()), os_receive(n_players);
    vector<gfpMatrix> shares_prng(n_players);

    // Request part of input.
    input_block_from_party_request_PRG(0, 0, first_blk_rows, os_send, shares_prng[0], os_receive[0]);
    for(int i = 1, startRow = first_blk_rows; i < n_players; i++, startRow += blk_rows){
        input_block_from_party_request_PRG(i, startRow, blk_rows, os_send, shares_prng[i], os_receive[i]);
    }

    // Wait part of input.
    finish_input_block_from_party_PRG(0, 0, first_blk_rows, os_send, shares_prng[0], os_receive[0]);
    for(size_t i = 1, startRow = first_blk_rows; i < n_players; i++, startRow += blk_rows){
        finish_input_block_from_party_PRG(i, startRow, blk_rows, os_send, shares_prng[i], os_receive[i]);
    }

    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}
/*********************************************************************************************
 * 
 *       Reveal Methods
 * 
 * There are serveral reveal methods with tiny difference.
 * The only difference is the range of the reveal.
 * 
 * 1. reveal_to_party(player_no): The specific party collects shares from the reconstuction set and calcuate the secret.
 * 2. reveal_to_Pking() : The specific party is Pking.
 * 3. reveal(): reveal_to_Ping at first, then Pking send the secret to other parties.
 * 
 * 
 * *******************************************************************************************/
void ShareBundle::reveal_to_party(int player_no)
{
    if(P->my_num()!=player_no && is_in_reconstruction_set(degree)){
        send_shares(player_no);
    }

    if(P->my_num()==player_no){
        reconstruct_secrets();
    }
}

void ShareBundle::reveal_to_pking()
{
    reveal_to_party(Pking);
}

gfpMatrix ShareBundle::reveal(bool active)
{
#ifdef DISPERSE_PKING
    // Dispersed version.
    return reveal_dispersed(active);
#else
    if(active && (degree == threshold)){  // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
        compare_view();
    }

    // Normal version.
    reveal_to_pking();
    octetStream o;
    if(P->my_num()==Pking){
        send_secrets(o);
    }else{
        receive_secrets(Pking, o);
    }

    // TODO_ZYS: the view of secrets need to be compared.
    if(active) {
        P->update_hash_state(o.get_data(), o.get_length());
    }

    if(active && (degree == threshold)){  // Only t-sharings need to be checked!
        RevealVerifier::push_matrix(secrets, shares);
    }

    return secrets;
#endif
}

// Since this reveal will check the correctness by Lagrange interpolation, we do not need to compare view.
void ShareBundle::reveal_to_party_with_check(int player_no)
{
    assert(degree == threshold);
    compare_view(); // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
    if(P->my_num()!=player_no){
        send_shares(player_no);
    }

    if(P->my_num()==player_no){
        // reconstruct_secret()
        octetStreams os(n_players);
        receive_shares_all(os);
    
        // calculate_secret()
        gfpMatrix sharings(n_players, shares.size());
        unpack_row(sharings, os);
        sharings.row(P->my_num()) = shares.reshaped<Eigen::RowMajor>().transpose();  // reshaped(): transfer to linear vector.
        // Determine whether n points are within a polynomial
        gfpMatrix reconstruct_points = sharings.block(0, 0, degree+1, sharings.cols());
        gfpMatrix checked_points = sharings.block(degree+1, 0, n_players-degree-1, sharings.cols()); // n-t-1 * shares.size()
        gfpMatrix computed_points = reconstruction_with_n_points_t * reconstruct_points;  // n-t * shares.size()
        secrets = computed_points.row(0).reshaped<Eigen::RowMajor>(secrets.rows(), secrets.cols());
        if (checked_points != computed_points.block(1, 0, n_players-degree-1, sharings.cols())) {
            throw reveal_error();
            // cout<<"ShareBundle ### reveal with check error: "<<endl;
            // assert(checked_points == computed_points.block(1, 0, n_players-degree-1, sharings.cols()));
        }
    }
}

gfpMatrix ShareBundle::reveal_with_check()
{
    assert(degree == threshold);
    compare_view(); // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
    // for(size_t player_no = 0; player_no < n_players; player_no++){
    //     if(P->my_num()!=player_no){
    //         send_shares(player_no);
    //     }
    // }

    // // reconstruct_secret()
    // octetStreams os(n_players);
    // for(size_t player_no = 0; player_no < n_players; player_no++)
    // if(P->my_num()!=player_no){
    //     // P->receive_player(player_no, os[player_no]);
    //     P->request_receive(player_no, os[player_no]); // Requst to receive from player_no
    // }

    // for(size_t player_no = 0; player_no < n_players; player_no++)
    // if(P->my_num()!=player_no){
    //     P->wait_receive(player_no, os[player_no]);
    // }

    send_shares_all();

    octetStreams os(n_players);
    receive_shares_all(os);

    // calculate_secret()
    gfpMatrix sharings(n_players, shares.size());
    unpack_row(sharings, os);
    sharings.row(P->my_num()) = shares.reshaped<Eigen::RowMajor>().transpose();  // reshaped(): transfer to linear vector.
    // Determine whether n points are within a polynomial
    gfpMatrix reconstruct_points = sharings.block(0, 0, degree+1, sharings.cols());
    gfpMatrix checked_points = sharings.block(degree+1, 0, n_players-degree-1, sharings.cols()); // n-t-1 * shares.size()
    gfpMatrix computed_points = reconstruction_with_n_points_t * reconstruct_points;  // n-t * shares.size()
    secrets = computed_points.row(0).reshaped<Eigen::RowMajor>(secrets.rows(), secrets.cols());
    if (checked_points != computed_points.block(1, 0, n_players-degree-1, sharings.cols())) {
        throw reveal_error();
        // cout<<"ShareBundle ### reveal with check error: "<<endl;
        // assert(checked_points == computed_points.block(1, 0, n_players-degree-1, sharings.cols()));
    }
    return secrets;
}

// Reveal to Pking and Pking truncate some bits, and send to other parties.
gfpMatrix ShareBundle::reveal_truncate(size_t precision, bool active)
{
#ifdef DISPERSE_PKING
    // Dispersed version.
    return reveal_truncate_dispersed(precision, active);
#else
    if(active && (degree == threshold)){  // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
        compare_view();
    }

    // Normal version
    reveal_to_pking();
    octetStream o;
    if(P->my_num()==Pking){
        send_secrets(o);
    }else{
        receive_secrets(Pking, o);
    }

    // TODO_ZYS: the view of secrets need to be compared.
    if(active) {
        P->update_hash_state(o.get_data(), o.get_length());
    }

    if(active && (degree == threshold)){  // Only t-sharings need to be checked!
        RevealVerifier::push_matrix(secrets, shares);
    }

    truncate_matrix(secrets, precision);

    return secrets;
#endif
}

// Different precision
gfpMatrix ShareBundle::reveal_truncate(vector<size_t> &precision, bool active)
{
#ifdef DISPERSE_PKING
    // Dispersed version.
    return reveal_truncate_dispersed(precision, active);
#else
    if(active && (degree == threshold)){  // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
        compare_view();
    }

    // Normal version
    reveal_to_pking();
    octetStream o;
    if(P->my_num()==Pking){
        send_secrets(o);
    }else{
        receive_secrets(Pking, o);
    }

    // TODO_ZYS: the view of secrets need to be compared.
    if(active) {
        P->update_hash_state(o.get_data(), o.get_length());
    }

    if(active && (degree == threshold)){  // Only t-sharings need to be checked!
        RevealVerifier::push_matrix(secrets, shares);
    }

    for(size_t i = 0; i < size(); i++){
        secrets(i).truncate(precision[i]);
    }

    return secrets;
#endif
}

// We disperse the Pking role to every party.
gfpMatrix ShareBundle::reveal_dispersed(bool active)
{
    if(active && (degree == threshold)){  // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
        compare_view();
    }

    reveal_blocks_dispersed();
    octetStream os_send;
    octetStreams os_recev(n_players);
    
    size_t startRow, nRows;
    partition_rows(startRow, nRows);
    // Send the corresponding block of secrets to other parties.
    pack_rows(secrets, startRow, nRows, os_send);
    P->request_send_all(os_send); //* BUG LOG  We cannot use P->send_all(os_send). Otherwise, each party will be waiting to send. Hence, the receivers need to request to receive after the senders request to send.
    
    P->request_receive_respective(os_recev);
    
    P->wait_send_all(os_send);
    // Receive other block of secrets from other parties.
    P->wait_receive_respective(os_recev);
    unpack_rows(secrets, os_recev);

    // TODO_ZYS: the view of secrets need to be compared.
    if(active) {
        os_recev[P->my_num()] = os_send;
        for(int i = 0; i < n_players; i++){
            P->update_hash_state(os_recev[i].get_data(), os_recev[i].get_length());
        }
        // The following may be used in large number of parties.
        // os_recev[P->my_num()] = os_send;
        // for(int i = 1; i < n_players; i++){
        //     os_recev[0].concat(os_recev[i]);
        // }
        // P->update_hash_state(os_recev[0].get_data(), os_recev[0].get_length());
    }

    if(active && (degree == threshold)){  // Only t-sharings need to be checked!
        RevealVerifier::push_matrix(secrets, shares);
    }

    return secrets;
}

// Compared to the above functionality, it truncates the secrets before sending to other paties.
gfpMatrix ShareBundle::reveal_truncate_dispersed(size_t precision, bool active)
{
    if(active && (degree == threshold)){  // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
        compare_view();
    }

    reveal_blocks_dispersed();
    size_t startRow, nRows;
    partition_rows(startRow, nRows);
    
    // Send the secrets to other parties.
    octetStream os_send;
    octetStreams os_recev(n_players);
    pack_rows(secrets, startRow, nRows, os_send);
    P->request_send_all(os_send);

    P->request_receive_respective(os_recev);
    
    P->wait_send_all(os_send);
    P->wait_receive_respective(os_recev);
    unpack_rows(secrets, os_recev);

    // TODO_ZYS: the view of secrets need to be compared.
    if(active) {
        os_recev[P->my_num()] = os_send;
        for(int i = 0; i < n_players; i++){
            P->update_hash_state(os_recev[i].get_data(), os_recev[i].get_length());
        }
    }

    if(active && (degree == threshold)){  // Only t-sharings need to be checked!
        RevealVerifier::push_matrix(secrets, shares);
    }

    // Truncate the secrets.
    truncate_matrix(secrets, precision);

    return secrets;
}

gfpMatrix ShareBundle::reveal_truncate_dispersed(vector<size_t> &precision, bool active)
{
    if(active && (degree == threshold)){  // Note_ZYS: When reveal t-sharings, must compare view first. Because the adversary may conduct differential attack.
        compare_view();
    }

    reveal_blocks_dispersed();
    size_t startRow, nRows;
    partition_rows(startRow, nRows);
    
    // Send the secrets to other parties.
    octetStream os_send;
    octetStreams os_recev(n_players);
    pack_rows(secrets, startRow, nRows, os_send);
    P->request_send_all(os_send);

    P->request_receive_respective(os_recev);
    
    P->wait_send_all(os_send);
    P->wait_receive_respective(os_recev);
    unpack_rows(secrets, os_recev);

    // TODO_ZYS: the view of secrets need to be compared.
    if(active) {
        os_recev[P->my_num()] = os_send;
        for(int i = 0; i < n_players; i++){
            P->update_hash_state(os_recev[i].get_data(), os_recev[i].get_length());
        }
    }

    if(active && (degree == threshold)){  // Only t-sharings need to be checked!
        RevealVerifier::push_matrix(secrets, shares);
    }

    // Truncate the secrets.
    for(size_t i = 0; i < size(); i++){
        secrets(i).truncate(precision[i]);
    }

    return secrets;
}

/******************************************************************************************
 * 
 *       Dispersed version of Input and Reveal
 * 
 * The motivation is that Pking is overloaded on reveal_to_pking and input_from_pking,
 * which always apperas in pairs in 'reduce_degree', 'truncate', and 'reduce_truncate'.
 * To ease the load of Pking, we dispersed the Pking role to every party.
 * 
 * We partition the matrix to '#players' blocks, each of which contains consecutive rows.
 * We implement the partition rule in 'partition_rows(startRow, nRows)'.
 * 
 * We implement two main functionality to disperse the operations of input and reveal.
 * One is input_blocks_dispersed(), the other is reveal_blocks_dispersed().
 * Each of them contains the request part and wait part.
 * 
 * - reveal_blocks_dispersed()
 *   - reveal_block_from_party_request()
 *   - finish_reveal_block_from_party()
 * - input_blocks_dispersed()
 *   - input_block_from_party_request()
 *   - finish_input_block_from_party()
 * 
 * ****************************************************************************************/

/**
 * @brief Partition rule: Partition the matrix to '#players' parts.
 * P0 takes charges of the first # rows.
 * P1 takes charges of the next # rows.
 * ...
 * The total number of rows could not be divided exactly by the '#players'.
 * We assign slightly more rows to P0.
 * 
 * @param startRow [out]
 * @param nRows [out]
 */
void ShareBundle::partition_rows(size_t &startRow, size_t &nRows)  // partition parameter: n_players
{
    size_t n_rows = rows() / n_players;
    size_t first_n_rows = rows() - n_rows * (n_players - 1);
    if(P->my_num()==0){
        startRow = 0;
        nRows = first_n_rows;
    }else{
        startRow = first_n_rows + n_rows * (P->my_num()-1);
        nRows = n_rows;
    }
}

/**
 * @brief There is no need to compare view and do reveal verifier in this function, as we do these checks outside this function.
 * 
 */
void ShareBundle::reveal_blocks_dispersed()
{
    // The reconstrunction only consists of degree+1 players, starting from P0.
    size_t n_relevant_players = degree + 1; 
    
    // Relevant players is the reconstruction set, containing P0, P1, ... P_{degree+1}.
    // They send corresponding block of shares to other parties.
    octetStreams os_send(n_players); // *BUG LOG: the size of octetStream should be n_players rather n_relevant_players.
    octetStreams os_receive(n_relevant_players); // *BUG LOG: the size of octetStream should be n_relevant_players.
    for(size_t i = 0; i < n_relevant_players; i++){
        reveal_block_from_party_request(i, os_send, os_receive[i]);  // os_receive is partition for every player
    }
    
    for(size_t i = 0; i < n_relevant_players; i++){
        finish_reveal_block_from_party(i, os_send, os_receive[i]);
    }

    // If P is in the reconstruction set, the shares of P are not sent to P because we avoid the loopback.
    // Hence, we need to fill our own shares in os_send into the os_receive.
    if(is_in_reconstruction_set(degree)){
        os_receive[P->my_num()] = os_send[P->my_num()];
    }
    
    // Each party calculate the corresponding block of secrets.
    size_t startRow, nRows;
    partition_rows(startRow, nRows);
    calculate_block_secrets(degree, startRow, nRows, os_receive);
    
}



/**
 * @brief Partition the matrix and reveal each block to corresponding party.
 * This is the request part of the above functionality.
 * Every party request to send blocks of shares to different party.
 * The receiver request to receive the corresponding block of shares.
 * It is the similar process with the input as every party distribute different block of shares to different party.
 * 
 * @param player_no Player ID who sends blocks of shares to different parties.
 * @param os_send Stores the different blocks of shares.
 * @param o_receive 
 */
void ShareBundle::reveal_block_from_party_request(int player_no, octetStreams &os_send, octetStream &o_receive)
{
    if(P->my_num() == player_no){
        pack_rows(shares, os_send);  // pack_rows is partition with parameter n_players
        P->request_send_respective(os_send);
    }else{P->request_receive(player_no, o_receive);}
}

// The wait part of the above functionality.
void ShareBundle::finish_reveal_block_from_party(int player_no, octetStreams &os_send, octetStream &o_receive)
{
    if(P->my_num() == player_no){
        P->wait_send_respective(os_send);
    }else{
        P->wait_receive(player_no, o_receive);
    }
}

/**
 * @brief Patition the matrix on row-grained.
 * Each party is responsible for one block.
 * The party sends the corresponding block of shares.
 * The receiver waits for the corresponding block of shares. 
 * 
 */
void ShareBundle::input_blocks_dispersed(bool active)
{
    size_t blk_rows = rows() / n_players;
    size_t first_blk_rows = rows() - blk_rows * (n_players - 1);
    octetStreams os_send(n_players), os_receive(n_players);

    // Request part of input.
    // P0 input d (=first_blk_rows) rows of the secrets, starting from row(0);
    input_block_from_party_request(0, 0, first_blk_rows, os_send, os_receive[0]);
    // Other parties input the next consecutive d'(=blk_rows) rows of the secrets.
    for(size_t i = 1, startRow = first_blk_rows; i < n_players; i++, startRow += blk_rows){
        input_block_from_party_request(i, startRow, blk_rows, os_send, os_receive[i]);
    }

    // Wait part of input.
    finish_input_block_from_party(0, 0, first_blk_rows, os_send, os_receive[0]);
    for(size_t i = 1, startRow = first_blk_rows; i < n_players; i++, startRow += blk_rows){
        finish_input_block_from_party(i, startRow, blk_rows, os_send, os_receive[i]);
    }

    if (active && (degree==threshold)) {
        // Correctness check
        InputVerifier::push_matrix(shares);
    }
}

/**
 * @brief Request part of the dispersed input.
 * 
 * @param player_no 
 * @param startRow start row of the block.
 * @param nRows number of rows in the block.
 * @param os_send octetStream to store the calculated sharings.
 * @param o_receive octetStream to receive the shares
 */
void ShareBundle::input_block_from_party_request(int player_no, const int &startRow, const int &nRows, octetStreams &os_send, octetStream &o_receive)
{
    if(player_no == P->my_num()){
        calculate_block_sharings(degree, startRow, nRows, os_send);
        P->request_send_respective(os_send);
    }else{P->request_receive(player_no, o_receive);}
}

// Wait part of the dispersed input.
void ShareBundle::finish_input_block_from_party(int player_no, const int &startRow, const int &nRows, octetStreams &os_send, octetStream &o_receive)
{
    if(player_no == P->my_num()){
        P->wait_send_respective(os_send);
    }else{
        P->wait_receive(player_no, o_receive);
        unpack_rows(shares, startRow, nRows, o_receive);
    }
}

/*********************************************************************
 * 
 *       Random Methods
 * 
 * *******************************************************************/
ShareBundle& ShareBundle::random()
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline()){
            Phase->generate_random_sharings(size());
        }
        else{
            Phase->switch_to_offline();
            Phase->generate_random_sharings(size());
            Phase->switch_to_online();
        }
    }
    RandomShare::get_randoms(RandomShare::queueRandom, shares);
    return *this;
}

/*********************************************************************
 * 
 *       Complicated Methods of ShareBundle
 *  ! Since the following op is complicated, we use the arg PhaseConfig to generate proper unbounded random sharings.
 *  ! If you can manage the unbounded mult random sharings, you can set trueOffline = true.
 * 
 * *******************************************************************/
/**
 * @brief Unbounded fan-in multiplication (prefix version).
 * For each row (a, b, c, d, ...)
 * Returns (a, ab, abc, abcd, ...)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_prefix_mult()
{
    DoubleShareBundle R(rows(), cols());
    R.unbounded_prefix_mult_random();

    ShareBundle res(rows(), cols());
    res.shares = shares.array() * R.aux_shares.array();  // Note_ZYS: shares can not be zero.

    // TODO_ZYS: Use PRZS
    // res.reduce_degree();
    
    // Here the PRZS is ommitted for simplicity. But it dosen't matter to the online phase.
    res.double_degree();
    res.reveal(true);

    // Verify
    MultVerifier::push_triples(shares, R.aux_shares, res.secrets);

    for(size_t i = 1; i < cols(); i++){
        res.secrets.col(i) = res.secrets.col(i).array() * res.secrets.col(i-1).array();
    }

    res.shares = res.secrets.array() * R.shares.array();
    return res;
}

ShareBundle ShareBundle::unbounded_postfix_mult()
{
    ShareBundle reverseShareBundle(rows(), cols());
    reverseShareBundle.shares = shares.rowwise().reverse();
    ShareBundle res = reverseShareBundle.unbounded_prefix_mult();
    res.shares.rowwise().reverseInPlace();
    return res;
}

/**
 * @brief Unbounded fan-in multiplication.
 * For each row (a, b, c, d, ...)
 * Returns abcd... (one entry)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_mult()
{
    DoubleShareBundle R(rows(), cols());
    R.unbounded_prefix_mult_random();
    

    ShareBundle res(rows(), cols());
    res.shares = shares.array() * R.aux_shares.array();
    // res.reduce_degree();// TODO: Use PRZS

    res.double_degree();
    res.reveal(true);

    // Verify
    MultVerifier::push_triples(shares, R.aux_shares, res.secrets);

    // * The following is different compared to the above functionality. 
    ShareBundle prod(rows(), 1);
    prod.shares = R.shares.rightCols(1);
    for(size_t i = 0; i < cols(); i++){
        prod.shares = prod.shares.array() * res.secrets.col(i).array();
    }

    return prod;
}

/**
 * @brief Evaluate the function value of each entry in the matrix.
 * 
 * @param phase If in true offline, we use the unbounded mult random sharings in the queue, and generate in separate offline otherwise.
 * @param degree 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::evalFunc(string fn, size_t degree)
{
    gfpVector funcConsts = getFuncConsts(fn, degree);
    ShareBundle expandMatrix(size(), degree); // Each row corresponds to one entry.
    for(size_t i = 0; i < size(); i++){
        expandMatrix.shares.row(i).setConstant(shares(i));
    }
    ShareBundle powMatrix = expandMatrix.unbounded_prefix_mult();// Calculate pow of each entry.

    gfpMatrix powMatrixPrime(size(), degree+1);
    powMatrixPrime.col(0).setConstant(1); // Pad 1 for the free coefficient.
    powMatrixPrime.middleCols(1, degree) = powMatrix.shares;
    
    ShareBundle res(rows(), cols());
    res.shares = (powMatrixPrime * funcConsts).reshaped<RowMajor>(rows(), cols()); // Evaluate the function given the coefficients.
    return res;
}


// Three layer multiplication: compute share(i.e., x) * y * z * w
// 1 round, 2 reveal
ShareBundle ShareBundle::three_layers_DN(const ShareBundle &y, const ShareBundle &z, const ShareBundle &w)
{
    assert(rows()==y.rows());
    assert(cols()==y.cols());
    assert(rows()==z.rows());
    assert(cols()==z.cols());
    assert(rows()==w.rows());
    assert(cols()==w.cols());

    DoubleShareBundlePair R_pairs(rows(), cols());
    R_pairs.mult_reduced_random_pair();

    ShareBundle combine(rows()<<1, cols());
    // Compute [u1]_2t = [xy]_2t + [r1]_2t and [u2]_2t = [zw]_2t + [r2]_2t
    combine.shares.topRows(rows()) = shares.array() * y.shares.array() + R_pairs.aux_shares.array();
    combine.shares.bottomRows(rows()) = z.shares.array() * w.shares.array() + R_pairs.other_aux_shares.array();
    combine.double_degree();
    combine.reveal(true);

    // Verify
    MultVerifier::push_triples(shares, y.shares, combine.secrets.topRows(rows()) - R_pairs.shares);
    MultVerifier::push_triples(z.shares, w.shares, combine.secrets.bottomRows(rows()) - R_pairs.other_shares);

    ShareBundle res(rows(), cols());
    res.shares = combine.secrets.topRows(rows()).array() * combine.secrets.bottomRows(rows()).array() // u1 * u2
                - combine.secrets.topRows(rows()).array() * R_pairs.other_shares.array() // u1 * [r2]_t
                - combine.secrets.bottomRows(rows()).array() * R_pairs.shares.array() // u2 * [r1]_t
                + R_pairs.prod_shares.array(); // [r1r2]_t
    return res;
}

// Reused three layer-mult: compute ([xyz1w1], [xyz2w2], ...) from input [x], [y] and ([z1], [z2], ...), ([w1], [w2], ...).
// Note that the size of vector is at most 4.
// 1 round, m+1 reveal.
vector<ShareBundle> ShareBundle::three_layers_DN(const ShareBundle &y, const vector<ShareBundle> &z, const vector<ShareBundle> &w)
{
    size_t len = z.size();
    assert(len==w.size());
    assert(0<len<=4);
    
    vector<ShareBundle> res_vector;
    if(len == 1) {
        res_vector.push_back(this->three_layers_DN(y, z[0], w[0]));
    } else {
        ReusedDoubleShareBundlePair R_pairs(rows(), cols());
        R_pairs.reused_mult_reduced_random_pairs();
        ShareBundle combine(rows()*(len+1), cols());
        // Compute [u1]_2t = [xy]_2t + [r1]_2t
        combine.shares.topRows(rows()) = shares.array() * y.shares.array() + R_pairs.aux_shares.array();
        for (int i = 0; i < len; i++)
        {
            ShareBundle tmp_z = z[i];
            ShareBundle tmp_w = w[i];
            assert(rows()==y.rows());
            assert(cols()==y.cols());
            assert(rows()==tmp_z.rows());
            assert(cols()==tmp_z.cols());
            assert(rows()==tmp_w.rows());
            assert(cols()==tmp_w.cols());
            
            // Compute [u2]_2t = [zw]_2t + [r2]_2t
            combine.shares(seqN((i+1)*rows(), rows()), seq(0, last)) = tmp_z.shares.array() * tmp_w.shares.array() + R_pairs.other_aux_shares(seqN(i*rows(), rows()), seq(0, last)).array();
        }
        combine.double_degree();
        combine.reveal(true);

        // Verify
        MultVerifier::push_triples(shares, y.shares, combine.secrets.topRows(rows()) - R_pairs.shares);
        for (int i = 0; i < len; i++)
        {
            MultVerifier::push_triples(z[i].shares, w[i].shares, combine.secrets(seqN((i+1)*rows(), rows()), seq(0, last)) - R_pairs.other_shares(seqN(i*rows(), rows()), seq(0, last)));
            ShareBundle res(rows(), cols());
            res.shares = combine.secrets.topRows(rows()).array() * combine.secrets(seqN((i+1)*rows(), rows()), seq(0, last)).array() // u1 * u2
                        - combine.secrets.topRows(rows()).array() * R_pairs.other_shares(seqN(i*rows(), rows()), seq(0, last)).array() // u1 * [r2]_t
                        - combine.secrets(seqN((i+1)*rows(), rows()), seq(0, last)).array() * R_pairs.shares.array() // u2 * [r1]_t
                        + R_pairs.prod_shares(seqN(i*rows(), rows()), seq(0, last)).array(); // [r1r2]_t
            
            res_vector.push_back(res);
        }  
    }
    return res_vector;
}

/**
 * @brief Unbounded fan-in multiplication with four elements by using 3L-DN protocol.
 * For each row (a, b, c, d)
 * Returns abcd (one entry)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_mult_4Elems() // 1 round, 2 reveal
{
    assert(shares.cols()==4);
    ShareBundle A(rows(), 1);
    A.shares = shares.col(0);
    ShareBundle B(rows(), 1);
    B.shares = shares.col(1);
    ShareBundle C(rows(), 1);
    C.shares = shares.col(2);
    ShareBundle D(rows(), 1);
    D.shares = shares.col(3);
    ShareBundle prod = A.three_layers_DN(B, C, D);

    return prod;
}

/**
 * @brief Reused unbounded fan-in multiplication with four elements by using 3L-DN protocol. Reuse the first 3 elements and at most compute 4 unbounded mult.
 * For each row (a, b, c, d1, d2, d3, d4)
 * Returns abcd1, abcd2, abcd3, abcd4
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_mult_4Elems_reused() // 1 round, m+1 reveal
{
    int num = shares.cols();
    assert(4<=num<=7);
    ShareBundle A(rows(), 1);
    A.shares = shares.col(0);
    ShareBundle B(rows(), 1);
    B.shares = shares.col(1);
    vector<ShareBundle> C;
    vector<ShareBundle> D;
    for(int i = 0; i < (num - 3); i++) {
        ShareBundle tmp_C(rows(), 1);
        tmp_C.shares = shares.col(2);
        C.push_back(tmp_C);
        ShareBundle tmp_D(rows(), 1);
        tmp_D.shares = shares.col(i+3);
        D.push_back(tmp_D);
    }
    vector<ShareBundle> res = A.three_layers_DN(B, C, D);
    ShareBundle prod(rows(), res.size());
    for(int i = 0; i < res.size(); i++) {
        prod.shares.col(i) = res[i].shares;
    }
    return prod;
}

/**
 * @brief Unbounded fan-in multiplication with three elements by using 2L-DN protocol.
 * For each row (a, b, c)
 * Returns abc (one entry)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_mult_3Elems() // 1 round, 2 reveal
{
    assert(shares.cols()==3);
    ShareBundle AB(rows(), 1);
    AB.shares = shares.col(0).array() * shares.col(1).array();
    ShareBundle C(rows(), 1);
    C.shares = shares.col(2);
    BeaverTriple triples(rows(), 1);
    AB.reduce_degree_1stLayer(C, triples);
    // Verify
    MultVerifier::push_triples(shares.col(0), shares.col(1), AB.shares);
    return triples.mult();
}

/**
 * @brief Reused unbounded fan-in multiplication with three elements by using 2L-DN protocol. Reuse the first 2 elements and at most compute 4 unbounded mult.
 * For each row (a, b, c1, c2, c3, c4)
 * Returns abc1, abc2, abc3, abc4
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_mult_3Elems_reused() // 1 round, m+1 reveal
{
    int num = shares.cols();
    assert(3<=num<=6);
    ShareBundle AB(rows(), 1);
    AB.shares = shares.col(0).array() * shares.col(1).array();
    vector<ShareBundle> C;
    vector<BeaverTriple> triples;
    for(int i = 0; i < (num - 2); i++) {
        ShareBundle tmp_C(rows(), 1);
        tmp_C.shares = shares.col(i+2);
        C.push_back(tmp_C);
        triples.push_back(BeaverTriple(rows(), 1));
    }
    AB.reduce_degree_1stLayer(C, triples);
    // Verify
    MultVerifier::push_triples(shares.col(0), shares.col(1), AB.shares);
    ShareBundle res(rows(), triples.size());
    for (int i = 0; i < triples.size(); i++)
    {
        res.shares.col(i) = triples[i].mult().shares;
    }
    
    return res;
}

/**
 * @brief Unbounded fan-in multiplication with four elements by using 3L-DN protocol (prefix version).
 * For each row (a, b, c, d)
 * Returns (a, ab, abc, abcd)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_prefix_mult_4Elems() // 1 round, 3 reveal
{
    assert(shares.cols()==4);
    DoubleShareBundlePair R_pairs(rows(), 1);
    R_pairs.mult_reduced_random_pair();
    DoubleShareBundle R(rows(), 1);
    R.reduced_random();  // [r_3]_t, [r_3]_2t

    ShareBundle combine(rows()*3, 1);
    // Compute [u1]_2t = [ab]_2t + [r1]_2t and [u2]_2t = [cd]_2t + [r2]_2t and [u3]_2t = [r1*(-c)]_2t + [r_3]_2t
    combine.shares.topRows(rows()) = shares.col(0).array() * shares.col(1).array() + R_pairs.aux_shares.array();
    combine.shares(seqN(rows(), rows()), seq(0, last)) = shares.col(2).array() * shares.col(3).array() + R_pairs.other_aux_shares.array();
    combine.shares.bottomRows(rows()) = (- shares.col(2)).array() * R_pairs.shares.array() + R.aux_shares.array();
    combine.double_degree();
    combine.reveal(true);

    // Verify
    MultVerifier::push_triples(shares.col(0), shares.col(1), combine.secrets.topRows(rows()) - R_pairs.shares);
    gfpMatrix verified_b(rows(), 2);
    verified_b << shares.col(3), R_pairs.shares;
    gfpMatrix verified_c(rows(), 2);
    verified_c << combine.secrets(seqN(rows(), rows()), seq(0, last)) - R_pairs.other_shares, R.shares - combine.secrets.bottomRows(rows());
    MultVerifier::push_reused_triples(shares.col(2), verified_b, verified_c);

    ShareBundle res(rows(), cols());
    res.shares.col(0) = shares.col(0);
    res.shares.col(1) = combine.secrets.topRows(rows()) - R_pairs.shares;
    res.shares.col(2) = combine.secrets.topRows(rows()).array() * shares.col(2).array() + combine.secrets.bottomRows(rows()).array() - R.shares.array();
    res.shares.col(3) = combine.secrets.topRows(rows()).array() * combine.secrets(seqN(rows(), rows()), seq(0, last)).array() // u1 * u2
                        - combine.secrets.topRows(rows()).array() * R_pairs.other_shares.array() // u1 * [r2]_t
                        - combine.secrets(seqN(rows(), rows()), seq(0, last)).array() * R_pairs.shares.array() // u2 * [r1]_t
                        + R_pairs.prod_shares.array(); // [r1r2]_t
    return res;
}

ShareBundle ShareBundle::unbounded_postfix_mult_4Elems() // 1 round, 3 reveal
{
    ShareBundle reverseShareBundle(rows(), cols());
    reverseShareBundle.shares = shares.rowwise().reverse();
    ShareBundle res = reverseShareBundle.unbounded_prefix_mult_4Elems();
    res.shares.rowwise().reverseInPlace();
    return res;
}

/**
 * @brief Unbounded fan-in multiplication with three elements by using 2L-DN protocol (prefix version).
 * For each row (a, b, c)
 * Returns (a, ab, abc)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_prefix_mult_3Elems() // 1 round, 2 reveal
{
    assert(shares.cols()==3);
    ShareBundle AB(rows(), 1);
    AB.shares = shares.col(0).array() * shares.col(1).array();
    ShareBundle C(rows(), 1);
    C.shares = shares.col(2);
    BeaverTriple triples(rows(), 1);
    AB.reduce_degree_1stLayer(C, triples);
    // Verify
    MultVerifier::push_triples(shares.col(0), shares.col(1), AB.shares);
    ShareBundle res(rows(), 3);
    res.shares.col(0) = shares.col(0);
    res.shares.col(1) = AB.shares;
    res.shares.col(2) = triples.mult().shares;
    return res;
}

ShareBundle ShareBundle::unbounded_postfix_mult_3Elems() // 1 round, 2 reveal
{
    ShareBundle reverseShareBundle(rows(), cols());
    reverseShareBundle.shares = shares.rowwise().reverse();
    ShareBundle res = reverseShareBundle.unbounded_prefix_mult_3Elems();
    res.shares.rowwise().reverseInPlace();
    return res;
}

/**
 * @brief Unbounded fan-in multiplication with mixed number of elements. The reused number is fixed to four.
 * For each row (a, b1, b2, b3, b4, 0, 0)
 * Returns (ab1, ab2, ab3, ab4)
 * For each row (a, b, c1, c2, c3, c4, 0)
 * Returns (abc1, abc2, abc3, abc4)
 * For each row (a, b, c, d1, d2, d3, d4)
 * Returns (abcd1, abcd2, abcd3, abcd4)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_mult_mixed_reused(size_t two_elem_rows, size_t three_elem_rows, size_t four_elem_rows)
{
    assert(shares.cols()==7);
    assert(shares.rows()==two_elem_rows+three_elem_rows+four_elem_rows);
    assert(two_elem_rows>0);
    assert(three_elem_rows>0);
    assert(four_elem_rows>0);

    // Compute two-elements mult
    DoubleShareBundle r(two_elem_rows, 4);
    r.reduced_random();
    gfpMatrix res1 = gfpMatrix(shares(seqN(0, two_elem_rows), 0).rowwise().replicate(4)).array() * shares(seqN(0, two_elem_rows), seqN(1, 4)).array() + r.aux_shares.array();

    // Compute three-elements mult
    DoubleShareBundle R(three_elem_rows, 5);
    R.reduced_random();
    gfpMatrix tmp_res2(three_elem_rows, 5);
    tmp_res2.col(0) = shares(seqN(two_elem_rows, three_elem_rows), 0).array() * shares(seqN(two_elem_rows, three_elem_rows), 1).array() + R.aux_shares.col(0).array();
    for (int i = 1; i < 5; i++)
    {
        tmp_res2.col(i) = (- shares(seqN(two_elem_rows, three_elem_rows), i+1)).array() * R.shares.col(0).array() + R.aux_shares.col(i).array();
    }

    // Compute four-elements mult
    ReusedDoubleShareBundlePair R_pairs(four_elem_rows, 1);
    R_pairs.reused_mult_reduced_random_pairs();
    gfpMatrix tmp_res3(four_elem_rows, 5);
    tmp_res3.col(0) = shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0).array()
                        * shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 1).array()
                        + R_pairs.aux_shares.array();
    for (int i = 1; i < 5; i++)
    {
        tmp_res3.col(i) = shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 2).array()
                            * shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), i+2).array()
                            + R_pairs.other_aux_shares(seqN((i-1)*four_elem_rows, four_elem_rows), 0).array();
    }

    // Batch open
    ShareBundle combine(rows(), 5);
    combine.shares.setConstant(0);
    combine.shares(seqN(0, two_elem_rows), seqN(0, 4)) = res1;
    combine.shares(seqN(two_elem_rows, three_elem_rows), seq(0, last)) = tmp_res2;
    combine.shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seq(0, last)) = tmp_res3;
    combine.double_degree();
    combine.reveal(true);

    // Set the result
    ShareBundle res(rows(), 4);
    // Set two-elements prods
    MultVerifier::push_reused_triples(shares(seqN(0, two_elem_rows), 0), 
                                shares(seqN(0, two_elem_rows), seqN(1, 4)),
                                combine.secrets(seqN(0, two_elem_rows), seqN(0, 4)) - r.shares);
    res.shares(seqN(0, two_elem_rows), seq(0, last)) = combine.secrets(seqN(0, two_elem_rows), seqN(0, 4)) - r.shares;
    // Set three-elements prods
    MultVerifier::push_triples(shares(seqN(two_elem_rows, three_elem_rows), 0),
                                shares(seqN(two_elem_rows, three_elem_rows), 1),
                                combine.secrets(seqN(two_elem_rows, three_elem_rows), 0) - R.shares.col(0));
    MultVerifier::push_reused_triples(R.shares.col(0), - shares(seqN(two_elem_rows, three_elem_rows), seqN(2, 4)),
                                combine.secrets(seqN(two_elem_rows, three_elem_rows), seqN(1, 4)) - R.shares(seqN(0, three_elem_rows), seqN(1, 4)));                     
    res.shares(seqN(two_elem_rows, three_elem_rows), seq(0, last)) =
                            gfpMatrix(combine.secrets(seqN(two_elem_rows, three_elem_rows), 0).rowwise().replicate(4)).array()
                            * shares(seqN(two_elem_rows, three_elem_rows), seqN(2, 4)).array() 
                            + combine.secrets(seqN(two_elem_rows, three_elem_rows), seqN(1, 4)).array()
                            - R.shares(seqN(0, three_elem_rows), seqN(1, 4)).array();
    // Set four-elements prods
    gfpMatrix R_pairs_other_shares(four_elem_rows, 4);
    R_pairs_other_shares.col(0) = R_pairs.other_shares(seqN(0, four_elem_rows), 0);
    R_pairs_other_shares.col(1) = R_pairs.other_shares(seqN(four_elem_rows, four_elem_rows), 0);
    R_pairs_other_shares.col(2) = R_pairs.other_shares(seqN(2*four_elem_rows, four_elem_rows), 0);
    R_pairs_other_shares.col(3) = R_pairs.other_shares(seqN(3*four_elem_rows, four_elem_rows), 0);
    gfpMatrix R_pairs_prod_shares(four_elem_rows, 4);
    R_pairs_prod_shares.col(0) = R_pairs.prod_shares(seqN(0, four_elem_rows), 0);
    R_pairs_prod_shares.col(1) = R_pairs.prod_shares(seqN(four_elem_rows, four_elem_rows), 0);
    R_pairs_prod_shares.col(2) = R_pairs.prod_shares(seqN(2*four_elem_rows, four_elem_rows), 0);
    R_pairs_prod_shares.col(3) = R_pairs.prod_shares(seqN(3*four_elem_rows, four_elem_rows), 0);
    MultVerifier::push_triples(shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0),
                                shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 1),
                                combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0) - R_pairs.shares);
    MultVerifier::push_reused_triples(shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 2),
                                shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(3, 4)),
                                combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(1, 4)) - R_pairs_other_shares);
    res.shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seq(0, last)) =
                            gfpMatrix(combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0).rowwise().replicate(4)).array()
                            * combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(1, 4)).array()
                            - gfpMatrix(combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0).rowwise().replicate(4)).array()
                            * R_pairs_other_shares.array()
                            - combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(1, 4)).array()
                            * gfpMatrix(R_pairs.shares.rowwise().replicate(4)).array()
                            + R_pairs_prod_shares.array();

    return res;
}

/**
 * @brief Unbounded fan-in multiplication with mixed number of elements. The reused number is fixed to sixteen.
 * For each row (a, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, 0, 0)
 * Returns (ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8, ab9, ab10, ab11, ab12, ab13, ab14, ab15, ab16)
 * For each row (a, b, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, 0)
 * Returns (abc1, abc2, abc3, abc4, abc5, abc6, abc7, abc8, abc9, abc10, abc11, abc12, abc13, abc14, abc15, abc16)
 * For each row (a, b, c, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16)
 * Returns (abcd1, abcd2, abcd3, abcd4, abcd5, abcd6, abcd7, abcd8, abcd9, abcd10, abcd11, abcd12, abcd13, abcd14, abcd15, abcd16)
 * 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::unbounded_mult_mixed_reused_16(size_t two_elem_rows, size_t three_elem_rows, size_t four_elem_rows)
{
    assert(shares.cols()==19);
    assert(shares.rows()==two_elem_rows+three_elem_rows+four_elem_rows);
    assert(two_elem_rows>0);
    assert(three_elem_rows>0);
    assert(four_elem_rows>0);

    // Compute two-elements mult
    DoubleShareBundle r_4_two(two_elem_rows, 16);
    r_4_two.reduced_random();
    gfpMatrix res1 = gfpMatrix(shares(seqN(0, two_elem_rows), 0).rowwise().replicate(16)).array() * shares(seqN(0, two_elem_rows), seqN(1, 16)).array() + r_4_two.aux_shares.array();

    // Compute three-elements mult
    DoubleShareBundle r_4_three(three_elem_rows, 17);
    r_4_three.reduced_random();
    gfpMatrix tmp_res2(three_elem_rows, 17);
    tmp_res2.col(0) = shares(seqN(two_elem_rows, three_elem_rows), 0).array() * shares(seqN(two_elem_rows, three_elem_rows), 1).array() + r_4_three.aux_shares.col(0).array();
    for (int i = 1; i < 17; i++)
    {
        tmp_res2.col(i) = (- shares(seqN(two_elem_rows, three_elem_rows), i+1)).array() * r_4_three.shares.col(0).array() + r_4_three.aux_shares.col(i).array();
    }

    // Compute four-elements mult
    ReusedDoubleShareBundlePair16 r_4_four(four_elem_rows, 1);
    r_4_four.reused_mult_reduced_random_pairs_16();
    gfpMatrix tmp_res3(four_elem_rows, 17);
    tmp_res3.col(0) = shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0).array()
                        * shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 1).array()
                        + r_4_four.aux_shares.array();
    for (int i = 1; i < 17; i++)
    {
        tmp_res3.col(i) = shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 2).array()
                            * shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), i+2).array()
                            + r_4_four.other_aux_shares(seqN((i-1)*four_elem_rows, four_elem_rows), 0).array();
    }

    // Batch open
    ShareBundle combine(rows(), 17);
    combine.shares.setConstant(0);
    combine.shares(seqN(0, two_elem_rows), seqN(0, 16)) = res1;
    combine.shares(seqN(two_elem_rows, three_elem_rows), seq(0, last)) = tmp_res2;
    combine.shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seq(0, last)) = tmp_res3;
    combine.double_degree();
    combine.reveal();

    // Set the result
    ShareBundle res(rows(), 16);
    // Set two-elements prods
    MultVerifier::push_reused_triples(shares(seqN(0, two_elem_rows), 0), 
                                shares(seqN(0, two_elem_rows), seqN(1, 16)),
                                combine.secrets(seqN(0, two_elem_rows), seqN(0, 16)) - r_4_two.shares);
    res.shares(seqN(0, two_elem_rows), seq(0, last)) = combine.secrets(seqN(0, two_elem_rows), seqN(0, 16)) - r_4_two.shares;
    // Set three-elements prods
    MultVerifier::push_triples(shares(seqN(two_elem_rows, three_elem_rows), 0),
                                shares(seqN(two_elem_rows, three_elem_rows), 1),
                                combine.secrets(seqN(two_elem_rows, three_elem_rows), 0) - r_4_three.shares.col(0));
    MultVerifier::push_reused_triples(r_4_three.shares.col(0), - shares(seqN(two_elem_rows, three_elem_rows), seqN(2, 16)),
                                combine.secrets(seqN(two_elem_rows, three_elem_rows), seqN(1, 16)) - r_4_three.shares(seqN(0, three_elem_rows), seqN(1, 16)));
    res.shares(seqN(two_elem_rows, three_elem_rows), seq(0, last)) =
                            gfpMatrix(combine.secrets(seqN(two_elem_rows, three_elem_rows), 0).rowwise().replicate(16)).array()
                            * shares(seqN(two_elem_rows, three_elem_rows), seqN(2, 16)).array() 
                            + combine.secrets(seqN(two_elem_rows, three_elem_rows), seqN(1, 16)).array()
                            - r_4_three.shares(seqN(0, three_elem_rows), seqN(1, 16)).array();
    // Set four-elements prods
    gfpMatrix r_4_four_other_shares(four_elem_rows, 16);
    gfpMatrix r_4_four_prod_shares(four_elem_rows, 16);
    for (size_t i = 0; i < 16; i++)
    {
        r_4_four_other_shares.col(i) = r_4_four.other_shares(seqN(i*four_elem_rows, four_elem_rows), 0);
        r_4_four_prod_shares.col(i) =  r_4_four.prod_shares(seqN(i*four_elem_rows, four_elem_rows), 0);
    }
    MultVerifier::push_triples(shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0),
                                shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 1),
                                combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0) - r_4_four.shares);
    MultVerifier::push_reused_triples(shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 2),
                                shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(3, 16)),
                                combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(1, 16)) - r_4_four_other_shares);
    res.shares(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seq(0, last)) =
                            gfpMatrix(combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0).rowwise().replicate(16)).array()
                            * combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(1, 16)).array()
                            - gfpMatrix(combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), 0).rowwise().replicate(16)).array()
                            * r_4_four_other_shares.array()
                            - combine.secrets(seqN(two_elem_rows+three_elem_rows, four_elem_rows), seqN(1, 16)).array()
                            * gfpMatrix(r_4_four.shares.rowwise().replicate(16)).array()
                            + r_4_four_prod_shares.array();
    return res;
}

ShareBundle ShareBundle::prefix_mult_log4_round()
{
    int blkNum = shares.cols() / 4;
    int batchSize = rows();
    #if defined(PR_31) 
        // Reshaped: four columns
        ShareBundle tmp_1((blkNum+1)*batchSize, 4);
        tmp_1.shares.setConstant(0);
        for (int i = 0; i < blkNum; i++)
        {
            tmp_1.shares(seqN(i*batchSize, batchSize), seq(0, last)) = shares(seq(0, last), seqN(4*i, 4));
        }
        tmp_1.shares(seqN(blkNum*batchSize, batchSize), seqN(0, 3)) = shares(seq(0, last), seqN(4*blkNum, 3));
        
        // First round: Prefix mult
        tmp_1 = tmp_1.unbounded_prefix_mult_4Elems();

        // Reshaped
        ShareBundle tmp_2(tmp_1.rows()-2*batchSize, 7);
        gfpMatrix tmp_3(batchSize, shares.cols());
        tmp_2.shares.setConstant(0);
        tmp_3.setConstant(0);
        // Set blk_1
        tmp_3(seq(0, last), seqN(0, 4)) = tmp_1.shares(seqN(0, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(16, 4)) = tmp_1.shares(seqN((batchSize<<2), batchSize), seq(0, last));
        // Set blk_2
        tmp_2.shares(seqN(0, batchSize), seqN(0, 5)) << tmp_1.shares.col(3)(seqN(0, batchSize)),
                                                        tmp_1.shares(seqN(batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(batchSize, batchSize), seqN(0, 5)) << tmp_1.shares.col(3)(seqN((batchSize<<2), batchSize)),
                                                                tmp_1.shares(seqN(5*batchSize, batchSize), seq(0, last));
        // Set blk_3
        tmp_2.shares(seqN(batchSize<<1, batchSize), seqN(0, 6)) << tmp_1.shares.col(3)(seqN(0, batchSize)),
                                                                tmp_1.shares.col(3)(seqN(batchSize, batchSize)),
                                                                tmp_1.shares(seqN((batchSize<<1), batchSize), seq(0, last));
        tmp_2.shares(seqN(3*batchSize, batchSize), seqN(0, 6)) << tmp_1.shares.col(3)(seqN(batchSize<<2, batchSize)),
                                                                tmp_1.shares.col(3)(seqN(5*batchSize, batchSize)),
                                                                tmp_1.shares(seqN(6*batchSize, batchSize), seq(0, last));
        // Set blk_4...
        tmp_2.shares(seqN(4*batchSize, batchSize), seq(0, last)) << tmp_1.shares.col(3)(seqN(0, batchSize)),
                                                                    tmp_1.shares.col(3)(seqN(batchSize, batchSize)),
                                                                    tmp_1.shares.col(3)(seqN(batchSize<<1, batchSize)),
                                                                    tmp_1.shares(seqN(3*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(5*batchSize, batchSize), seqN(0, 6)) << tmp_1.shares.col(3)(seqN(batchSize<<2, batchSize)),
                                                                    tmp_1.shares.col(3)(seqN(5*batchSize, batchSize)),
                                                                    tmp_1.shares.col(3)(seqN(6*batchSize, batchSize)),
                                                                    tmp_1.shares(seqN(7*batchSize, batchSize), seqN(0, 3));
        
        // Second round: Reused unbounded mult with mixed elements number
        tmp_2 = tmp_2.unbounded_mult_mixed_reused(batchSize<<1, batchSize<<1, batchSize<<1);
        // Resized into tmp_3
        // Set blk_2
        tmp_3(seq(0, last), seqN(4, 4)) = tmp_2.shares(seqN(0, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(20, 4)) = tmp_2.shares(seqN(batchSize, batchSize), seq(0, last));
        // Set blk_3
        tmp_3(seq(0, last), seqN(8, 4)) = tmp_2.shares(seqN((batchSize<<1), batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(24, 4)) = tmp_2.shares(seqN(3*batchSize, batchSize), seq(0, last));
        // Set blk_4
        tmp_3(seq(0, last), seqN(12, 4)) = tmp_2.shares(seqN(4*batchSize, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(28, 3)) = tmp_2.shares(seqN(5*batchSize, batchSize), seqN(0, 3));
        
        // Third round: Reused unbounded mult with two elements number
        ShareBundle res_tmp(batchSize, 15);
        res_tmp.shares = tmp_3.col(15).rowwise().replicate(15).array() * tmp_3(seq(0, last), seqN(16, 15)).array();
        res_tmp.reduce_degree();
        MultVerifier::push_reused_triples(tmp_3.col(15), tmp_3(seq(0, last), seqN(16, 15)), res_tmp.shares);

        // Reshaped
        ShareBundle res(batchSize, shares.cols());
        res.shares(seq(0, last), seqN(0, 16)) = tmp_3(seq(0, last), seqN(0, 16));
        res.shares(seq(0, last), seqN(16, 15)) = res_tmp.shares;

        return res;
         
    #elif defined(PR_61)
        // Reshaped: four columns
        gfpMatrix tmp_1((blkNum+1)*batchSize, 4);
        ShareBundle tmp_tmp_1(blkNum*batchSize, 4);
        tmp_1.setConstant(0);
        tmp_tmp_1.shares.setConstant(0);
        for (int i = 0; i < blkNum; i++)
        {
            tmp_tmp_1.shares(seqN(i*batchSize, batchSize), seq(0, last)) = shares(seq(0, last), seqN((i<<2), 4));
        }
        // First round: Prefix mult
        tmp_1(seqN(0, blkNum*batchSize), seq(0, last)) = tmp_tmp_1.unbounded_prefix_mult_4Elems().shares;
        tmp_1.col(0)(seqN(blkNum*batchSize, batchSize)) = shares.col(60)(seq(0, last));

        // Reshaped
        ShareBundle tmp_2(tmp_1.rows()-(batchSize<<2), 7);  // From blk_2...
        gfpMatrix tmp_3(batchSize, shares.cols());
        tmp_2.shares.setConstant(0);
        tmp_3.setConstant(0);
        // Set blk_1
        tmp_3(seq(0, last), seqN(0, 4)) = tmp_1(seqN(0, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(16, 4)) = tmp_1(seqN((batchSize<<2), batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(32, 4)) = tmp_1(seqN((batchSize<<3), batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(48, 4)) = tmp_1(seqN((batchSize*12), batchSize), seq(0, last));
        // Set blk_2
        tmp_2.shares(seqN(0, batchSize), seqN(0, 5)) << tmp_1.col(3)(seqN(0, batchSize)),
                                                        tmp_1(seqN(batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(batchSize, batchSize), seqN(0, 5)) << tmp_1.col(3)(seqN((batchSize<<2), batchSize)),
                                                                tmp_1(seqN(5*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN((batchSize<<1), batchSize), seqN(0, 5)) << tmp_1.col(3)(seqN((batchSize<<3), batchSize)),
                                                                tmp_1(seqN(9*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(3*batchSize, batchSize), seqN(0, 5)) << tmp_1.col(3)(seqN((batchSize*12), batchSize)),
                                                                tmp_1(seqN(batchSize*13, batchSize), seq(0, last));
        // Set blk_3
        tmp_2.shares(seqN((batchSize<<2), batchSize), seqN(0, 6)) << tmp_1.col(3)(seqN(0, batchSize)),
                                                                tmp_1.col(3)(seqN(batchSize, batchSize)),
                                                                tmp_1(seqN((batchSize<<1), batchSize), seq(0, last));
        tmp_2.shares(seqN(5*batchSize, batchSize), seqN(0, 6)) << tmp_1.col(3)(seqN(batchSize<<2, batchSize)),
                                                                tmp_1.col(3)(seqN(5*batchSize, batchSize)),
                                                                tmp_1(seqN(6*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(6*batchSize, batchSize), seqN(0, 6)) << tmp_1.col(3)(seqN((batchSize<<3), batchSize)),
                                                                tmp_1.col(3)(seqN(9*batchSize, batchSize)),
                                                                tmp_1(seqN(10*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(7*batchSize, batchSize), seqN(0, 6)) << tmp_1.col(3)(seqN((batchSize*12), batchSize)),
                                                                tmp_1.col(3)(seqN(13*batchSize, batchSize)),
                                                                tmp_1(seqN(14*batchSize, batchSize), seq(0, last));
        // Set blk_4...
        tmp_2.shares(seqN((batchSize<<3), batchSize), seq(0, last)) << tmp_1.col(3)(seqN(0, batchSize)),
                                                                    tmp_1.col(3)(seqN(batchSize, batchSize)),
                                                                    tmp_1.col(3)(seqN(batchSize<<1, batchSize)),
                                                                    tmp_1(seqN(3*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(9*batchSize, batchSize), seq(0, last)) << tmp_1.col(3)(seqN(batchSize<<2, batchSize)),
                                                                    tmp_1.col(3)(seqN(5*batchSize, batchSize)),
                                                                    tmp_1.col(3)(seqN(6*batchSize, batchSize)),
                                                                    tmp_1(seqN(7*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(10*batchSize, batchSize), seq(0, last)) << tmp_1.col(3)(seqN(batchSize<<3, batchSize)),
                                                                    tmp_1.col(3)(seqN(9*batchSize, batchSize)),
                                                                    tmp_1.col(3)(seqN(10*batchSize, batchSize)),
                                                                    tmp_1(seqN(11*batchSize, batchSize), seq(0, last));
        tmp_2.shares(seqN(11*batchSize, batchSize), seq(0, last)) << tmp_1.col(3)(seqN(12*batchSize, batchSize)),
                                                                    tmp_1.col(3)(seqN(13*batchSize, batchSize)),
                                                                    tmp_1.col(3)(seqN(14*batchSize, batchSize)),
                                                                    tmp_1(seqN(15*batchSize, batchSize), seq(0, last));
        
        // Second round: Reused unbounded mult with mixed elements number
        tmp_2 = tmp_2.unbounded_mult_mixed_reused(batchSize<<2, batchSize<<2, batchSize<<2);
        // Resized into tmp_3
        // Set blk_2
        tmp_3(seq(0, last), seqN(4, 4)) = tmp_2.shares(seqN(0, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(20, 4)) = tmp_2.shares(seqN(batchSize, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(36, 4)) = tmp_2.shares(seqN((batchSize<<1), batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(52, 4)) = tmp_2.shares(seqN(3*batchSize, batchSize), seq(0, last));
        // Set blk_3
        tmp_3(seq(0, last), seqN(8, 4)) = tmp_2.shares(seqN((batchSize<<2), batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(24, 4)) = tmp_2.shares(seqN(5*batchSize, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(40, 4)) = tmp_2.shares(seqN(6*batchSize, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(56, 4)) = tmp_2.shares(seqN(7*batchSize, batchSize), seq(0, last));
        // Set blk_4
        tmp_3(seq(0, last), seqN(12, 4)) = tmp_2.shares(seqN((batchSize<<3), batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(28, 4)) = tmp_2.shares(seqN(9*batchSize, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(44, 4)) = tmp_2.shares(seqN(10*batchSize, batchSize), seq(0, last));
        tmp_3(seq(0, last), seqN(60, 1)) = tmp_2.shares(seqN(11*batchSize, batchSize), seqN(0, 1));
        
        // Third round: Reused unbounded mult mixed elements number
        ShareBundle res_tmp(batchSize*3, 19);
        res_tmp.shares.setConstant(0);
        // Set reused unbounded mult with two elements number
        res_tmp.shares(seqN(0, batchSize), seqN(0, 17)) << tmp_3.col(15),
                                                        tmp_3(seq(0, last), seqN(16, 16));
        // Set reused unbounded mult with three elements number
        res_tmp.shares(seqN(batchSize, batchSize), seqN(0, 18)) << tmp_3.col(15),
                                                                tmp_3.col(31),
                                                                tmp_3(seq(0, last), seqN(32, 16));
        // Set reused unbounded mult with four elements number
        res_tmp.shares(seqN((batchSize<<1), batchSize), seqN(0, 16)) << tmp_3.col(15),
                                                                        tmp_3.col(31),
                                                                        tmp_3.col(47),
                                                                        tmp_3(seq(0, last), seqN(48, 13));
        res_tmp = res_tmp.unbounded_mult_mixed_reused_16(batchSize, batchSize, batchSize);

        // Reshaped
        ShareBundle res(batchSize, shares.cols());
        res.shares(seq(0, last), seqN(0, 16)) = tmp_3(seq(0, last), seqN(0, 16));
        res.shares(seq(0, last), seqN(16, 16)) = res_tmp.shares(seqN(0, batchSize), seq(0, last));
        res.shares(seq(0, last), seqN(32, 16)) = res_tmp.shares(seqN(batchSize, batchSize), seq(0, last));
        res.shares(seq(0, last), seqN(48, 13)) = res_tmp.shares(seqN((batchSize<<1), batchSize), seqN(0, 13));

        return res;
        
    #endif
}

ShareBundle ShareBundle::postfix_mult_log4_round()
{
    ShareBundle reverseShareBundle(rows(), cols());
    reverseShareBundle.shares = shares.rowwise().reverse();
    ShareBundle res = reverseShareBundle.prefix_mult_log4_round();
    res.shares.rowwise().reverseInPlace();
    return res;
}


// Get the least significant bit. (use an edaBit)
BitBundle ShareBundle::get_LSB()const
{
    BitBundle rBits(size());// size() by BITS_LENGTH
    ShareBundle rField(rows(), cols());
    rBits.solved_random(rField);

    BitBundle rLSB(rows(), cols());
    rLSB.shares = rBits.shares.col(0).reshaped<RowMajor>(rows(), cols());

    ShareBundle mask(rows(), cols());
    mask.shares = shares + rField.shares;
    mask.reveal(true);

    gfpMatrix maskLSB(rows(), cols());
    for(size_t i = 0; i < size(); i++){
        maskLSB(i) = mask.secret()(i) & 1;
    }

    BitBundle lsb = bitwise_xor(rLSB, maskLSB);

    mask.resize(rows()*cols(), 1);// ->vector.

    gfpMatrix maskBits(rows()*cols(), BITS_LENGTH);  // The size of Bitwise matrix is rows()*cols()=size(), BITS_LENGTH
    
    // decompose_bits(mask.secret(), BITS_LENGTH, maskBits);
    
    // First method
    // BitBundle is_wrap = less_than_unsigned(maskBits, rBits);

    // Second method
    // The complexity is 5 round when b is known in (a<b)
    // BUG LOG But When a = 0, MSB(a) = 1
    // BitBundle is_wrap = less_than_unsigned(rBits, maskBits);
    // is_wrap.shares = 1 - is_wrap.shares.array();

    // A Third way: compute LT(p-rBits, p-masked)
    decompose_bits(-mask.secret(), BITS_LENGTH, maskBits);  // When using Mersenne prime, -a = p-a, -a is equal to the bitwise inversion of a
    rBits.shares = 1 - rBits.shares.array();
    BitBundle is_wrap = less_than_unsigned(rBits, maskBits);

    is_wrap.resize(rows(), cols());
    
    return bitwise_xor(lsb, is_wrap);
}

BitBundle ShareBundle::get_MSB()const
{
    ShareBundle res(rows(), cols());
    res.shares = 2 * shares.array();
    return res.get_LSB();
}

// naive method
BitBundle ShareBundle::deltaReLU()const
{
    BitBundle res(rows(), cols());
    res.shares.setConstant(1);
    res.shares -= get_MSB().shares;
    return res;
}

// naive method
ShareBundle ShareBundle::ReLU()const
{
    ShareBundle res(rows(), cols());
    gfpMatrix tmp = deltaReLU().shares;
    res.shares = shares.array() * tmp.array();
    res.reduce_degree();

    // Verify
    MultVerifier::push_triples(shares, tmp, res.shares);

    return res;
}

// Get ReLUPrime and ReLU (naive method)
void ShareBundle::ReLU(ShareBundle &reluPrime, ShareBundle &relu)const
{
    reluPrime = deltaReLU();
    relu.shares = shares.array() * reluPrime.shares.array();
    relu.reduce_degree();

    // Verify
    MultVerifier::push_triples(shares, reluPrime.shares, relu.shares);

    return;
}

// * Test for two rounds ReLU.
// More elegant implementation is referred to ReLU_opt
void ShareBundle::ReLU_opt_test(ShareBundle &deltaReLU, ShareBundle &relu)const
{
    ShareBundle x(rows(), cols());
    x.shares = 2 * shares.array();

    // x0 (LSB) = x0_prime * x0_prime
    ShareBundle x0_prime = x.get_LSB_impared();

    // Two layer multiplication:
    // 1. The first layer: compute the last layer of mult in LSB circuit
    // 2. The second layer: only the first input wire is from the output of the first layer.
    ShareBundle x0(rows(), cols());
    x0.shares = x0_prime.shares.array() * x0_prime.shares.array();

    // This triple is to compute [x0]_t * [y]_t where [x0]_t = u - [a]_t
    BeaverTriple triple(rows(), cols());
    x0.reduce_degree_1stLayer(*this, triple);
    MultVerifier::push_triples(x0_prime.shares, x0_prime.shares, x0.shares);

    // deltaReLU = 1 - MSB(x') = 1 - LSB(x) where x = 2x' and MSB(x')=LSB(2x')
    // What we actually compute is (1-[x0]_t)*[y]_t where (1-[x0]_t) = (1-u) - (-[a]_t)

    // 2. The second layer:
    // The first input wire changes from [x0] to (-[x0]) + 1.
    deltaReLU.shares = 1 - x0.shares.array();
    triple.x_times(-gfpScalar(1)).x_plus(1); // (-x) + 1
    
    relu = triple.mult();
    return;
}

// Get the least significant bit without the last xor step.
// We remove the last layer of multiplication in LSB circuit.
// We want to compose a two-layer multiplication, which can be computed in one round, rather than two rounds.
ShareBundle ShareBundle::get_LSB_impared()const
{
    BitBundle rBits(size());// size() by BITS_LENGTH
    ShareBundle rField(rows(), cols());
    rBits.solved_random(rField);

    BitBundle rLSB(rows(), cols());
    rLSB.shares = rBits.shares.col(0).reshaped<RowMajor>(rows(), cols());

    ShareBundle mask(rows(), cols());
    mask.shares = shares + rField.shares;
    mask.reveal(true);

    gfpMatrix maskLSB(rows(), cols());
    for(size_t i = 0; i < size(); i++){
        maskLSB(i) = mask.secret()(i) & 1;
    }

    BitBundle lsb = bitwise_xor(rLSB, maskLSB);  // [b]_p in paper

    mask.resize(rows()*cols(), 1);// ->vector.
    gfpMatrix maskBits(rows()*cols(), BITS_LENGTH);
    
    decompose_bits(-mask.secret(), BITS_LENGTH, maskBits);
    rBits.shares = 1 - rBits.shares.array();

    BitBundle is_wrap = less_than_unsigned(rBits, maskBits);  // [c]_p in paper


    is_wrap.resize(rows(), cols());
    
    // The original return.
    // return bitwise_xor(lsb, is_wrap);

    // New return which eliminates the last multiplication.
    ShareBundle ret_diff(rows(), cols());
    ret_diff.shares = lsb.shares - is_wrap.shares;
    return ret_diff;  // [b-c]_p in paper
}

BitBundle ShareBundle::get_LSB_opt(vector<ShareBundle> &y, vector<BeaverTriple> &triples)const
{
    ShareBundle x0_prime = get_LSB_impared();

    // 1. The first layer of mult.
    BitBundle x0(rows(), cols());
    x0.shares = x0_prime.shares.array() * x0_prime.shares.array(); // reduce the degree
    x0.reduce_degree_1stLayer(y, triples); // compute the BeaverTriple for the second layer of mult

    // Verify
    MultVerifier::push_triples(x0_prime.shares, x0_prime.shares, x0.shares);
    // funcCheckMaliciousWiseMul(x0_prime.shares, x0_prime.shares, x0.shares);

    return x0; // [(b-c)*(b-c)]_t = b \xor c = LSB(a) in paper.
}

BitBundle ShareBundle::get_MSB_opt(vector<ShareBundle> &y, vector<BeaverTriple> &triples)const
{
    // MSB(x) = LSB(2x)
    ShareBundle res(rows(), cols());
    res.shares = 2 * shares.array();
    return res.get_LSB_opt(y, triples);
}

BitBundle ShareBundle::deltaReLU_opt(vector<ShareBundle> &y, vector<BeaverTriple> &triples)const
{
    BitBundle res = get_MSB_opt(y, triples);  // [(b-c)*(b-c)]_t = b \xor c in paper.
    // res = 1 - res
    res.shares = 1 - res.shares.array();  // DReLU(a) = 1 - b \xor c in paper.
    for(size_t i = 0; i < triples.size(); i++){
        triples[i].x_times(-gfpScalar(1)).x_plus(1);  // Do the same operations in triples.
    }
    return res;
}

// Optimize ReLU with 2-layer DN multiplication;
void ShareBundle::ReLU_opt(ShareBundle &deltaReLU, ShareBundle &relu)const
{
    vector<ShareBundle> y{*this};

    vector<BeaverTriple> triples{BeaverTriple(rows(), cols())};
    
    deltaReLU.shares = deltaReLU_opt(y, triples).shares;
    relu.shares = triples[0].mult().shares;
    return;
}

// Only return ReLU
ShareBundle ShareBundle::ReLU_opt()const
{
    ShareBundle reluPrime(rows(), cols()), relu(rows(), cols());
    ReLU_opt(reluPrime, relu);
    return relu;
}

// * Max function

// The sequential method to get MaxPool and MaxPoolPrime.
// Time is propotional to length of max.
void ShareBundle::seqMaxpoolRowwise(ShareBundle &maxRes, ShareBundle &maxIdx)const
{
    assert(maxRes.cols()==1);
    // maxRes.shares = shares.col(0);
    // maxIdx.shares.setConstant(0);
    // maxIdx.shares.col(0).setConstant(1);

    ShareBundle maxTmp(rows(), cols()+1);
    maxTmp.shares.col(0) = shares.col(0); //It stores the max value of each row.
    maxTmp.shares.rightCols(cols()).setConstant(0);// It stores the max idx vector.
    maxTmp.shares.col(1).setConstant(1);

    // For verify
    gfpMatrix aTmp(rows(), cols()+1);
    gfpMatrix bTmp(rows(), cols()+1);
    gfpMatrix cTmp(rows(), cols()+1);

    for(size_t i = 1; i < cols(); i++){
        ShareBundle deltaMax(rows(), 1);
        deltaMax.shares = maxTmp.shares.col(0) - shares.col(i);
        
        gfpMatrix currentIdx(1, cols());
        currentIdx.setConstant(0);
        currentIdx(0, i) = 1;

        ShareBundle cond = deltaMax.deltaReLU(); // cond = 0 : need to update
        aTmp.col(0) = cond.shares;
        bTmp.col(0) =  maxTmp.shares.col(0) - shares.col(i);
        // Update the max value
        maxTmp.shares.col(0) = (1 - cond.shares.array()) * shares.col(i).array() + cond.shares.array() * maxTmp.shares.col(0).array();

        // Update the max idx.
        for(size_t j = 0; j < rows(); j++){
            aTmp.row(j).tail(cols()).setConstant(cond.shares(j, 0));
            bTmp.row(j).tail(cols()) = maxTmp.shares.row(j).tail(cols());
            maxTmp.shares.row(j).tail(cols()) = (1 - cond.shares(j, 0)) * currentIdx.array() + cond.shares(j, 0) * maxTmp.shares.row(j).tail(cols()).array();
        }

        maxTmp.reduce_degree();
        cTmp.col(0) = maxTmp.shares.col(0) - shares.col(i);
        for(size_t j = 0; j < rows(); j++){
            cTmp.row(j).tail(cols()) = maxTmp.shares.row(j).tail(cols()).array() - ((1 - cond.shares(j, 0)) * currentIdx.array());
        }
        MultVerifier::push_triples(aTmp, bTmp, cTmp);
    }

    maxRes.shares = maxTmp.shares.col(0);
    maxIdx.shares = maxTmp.shares.rightCols(cols());
}

/**
 * @brief Use the hierachical way to get the max value and its corresponding index.
 * Compared to the above implemention, this only need log(n) round max(,).
 * @param depth starts from 1
 * @param nCols 
 * @param origin [in] The original share bundles.
 * @param maxIdx [out] The one-hot index. Note: the maxIdx must initialized with all ones vector.
 * @return ShareBundle [out] The max value of each row.
 */
ShareBundle ShareBundle::hierMaxpoolRowwise(int depth, size_t nCols, const ShareBundle &origin, ShareBundle &maxIdx)const
{
    size_t ySize = nCols/2;
    if(nCols%2==1){ySize++;}
    int blkSize = 1<<depth;// depth starts from 1

    ShareBundle blkBundle(rows(), nCols/2);
    for(size_t i = 0; i < nCols/2; i++){
        blkBundle.shares.col(i) = shares.col(i*2) - shares.col(i*2+1);
    }
    ShareBundle cond = blkBundle.deltaReLU(); // cond = 1, left > right -> max = (i*2)
    // cout<<cond.reveal()<<endl;

    ShareBundle maxBundle(rows(), ySize);
    //For verify
    gfpMatrix aTmp(rows(), nCols/2);
    gfpMatrix bTmp(rows(), nCols/2);
    gfpMatrix cTmp(rows(), nCols/2);
    int id_num = ((nCols/2) * blkSize) > origin.cols() ? origin.cols() : ((nCols/2) * blkSize);
    gfpMatrix aTmpIdx(rows(), id_num);
    gfpMatrix bTmpIdx(rows(), id_num);
    gfpMatrix cTmpIdx(rows(), id_num);
    for(size_t i = 0; i < nCols/2; i++){
        // max(a_{i*2}, a_{i*2+1}) = cond_i * a_{i*2} + (1-cond_i) * a_{i*2+1}
        aTmp.col(i) = cond.shares.col(i);
        bTmp.col(i) = shares.col(i*2) - shares.col(i*2+1);
        maxBundle.shares.col(i) = cond.shares.col(i).array() * shares.col(i*2).array() + (1 - cond.shares.col(i).array()) * shares.col(i*2+1).array();
        int start = i*blkSize;
        int leftBlkSize = blkSize/2;
        int rightBlkSize = (start+blkSize) > origin.cols() ? (origin.cols()-start-leftBlkSize) : leftBlkSize;
        // cout<<leftBlkSize<<" "<<rightBlkSize<<endl;
        // gfpMatrix(col.rowwise().replicate(N)): replicate N cols
        // maxIdx.shares.middleCols(start, leftBlkSize) = cond.shares.col(i).array() * maxIdx.shares.middleCols(start, leftBlkSize).array();
        aTmpIdx.middleCols(start, leftBlkSize) = gfpMatrix(cond.shares.col(i).rowwise().replicate(leftBlkSize));
        bTmpIdx.middleCols(start, leftBlkSize) = maxIdx.shares.middleCols(start, leftBlkSize);
        maxIdx.shares.middleCols(start, leftBlkSize) = gfpMatrix(cond.shares.col(i).rowwise().replicate(leftBlkSize)).array() * maxIdx.shares.middleCols(start, leftBlkSize).array();
        // maxIdx.shares.middleCols(start+leftBlkSize, rightBlkSize) = (1 - cond.shares.col(i).array()) * maxIdx.shares.middleCols(start+leftBlkSize, rightBlkSize).array();
        aTmpIdx.middleCols(start+leftBlkSize, rightBlkSize) = gfpMatrix((1 - cond.shares.col(i).array()).rowwise().replicate(rightBlkSize));
        bTmpIdx.middleCols(start+leftBlkSize, rightBlkSize) = maxIdx.shares.middleCols(start+leftBlkSize, rightBlkSize);
        maxIdx.shares.middleCols(start+leftBlkSize, rightBlkSize) = gfpMatrix((1 - cond.shares.col(i).array()).rowwise().replicate(rightBlkSize)).array() * maxIdx.shares.middleCols(start+leftBlkSize, rightBlkSize).array();
    }
    // maxBundle.reduce_degree();
    // maxIdx.reduce_degree();

    // Reduce degree at the same time.
    ShareBundle tmpMatrix(rows(), maxBundle.cols()+maxIdx.cols());
    tmpMatrix.shares<<maxBundle.shares, maxIdx.shares;
    tmpMatrix.reduce_degree();
    maxBundle.shares = tmpMatrix.shares.leftCols(maxBundle.cols());
    maxIdx.shares = tmpMatrix.shares.rightCols(maxIdx.cols());
    // Verify
    for(size_t i = 0; i < nCols/2; i++){
        cTmp.col(i) = maxBundle.shares.col(i) - shares.col(i*2+1);
        int start = i*blkSize;
        int leftBlkSize = blkSize/2;
        int rightBlkSize = (start+blkSize) > origin.cols() ? (origin.cols()-start-leftBlkSize) : leftBlkSize;
        cTmpIdx.middleCols(start, leftBlkSize) = maxIdx.shares.middleCols(start, leftBlkSize);
        cTmpIdx.middleCols(start+leftBlkSize, rightBlkSize) = maxIdx.shares.middleCols(start+leftBlkSize, rightBlkSize);
    }
    MultVerifier::push_triples(aTmp, bTmp, cTmp);
    MultVerifier::push_triples(aTmpIdx, bTmpIdx, cTmpIdx);
    aTmp.resize(0, 0);
    bTmp.resize(0, 0);
    cTmp.resize(0, 0);
    aTmpIdx.resize(0, 0);
    bTmpIdx.resize(0, 0);
    cTmpIdx.resize(0, 0);

    if(nCols%2==1){
        maxBundle.shares.rightCols(1) = shares.rightCols(1);
    }

    // Stop condition: Get the maximum.
    if(ySize==1){
        return maxBundle;
    }
    return maxBundle.hierMaxpoolRowwise(depth+1, ySize, origin, maxIdx);
}

/**
 * @brief Use the two-layer multiplication to optmize the round complexity by factor 2.
 * The functionality is the same as 'hierMaxpoolRowwise'.
 * But we will not evaluate the 'cond' = deltaReLU immediately.
 * Instead we will prepare all the Beaver triples of the second layer-multiplication, so
 * the input for the second layer should handled before the real evaluation of deltaReLU.
 * 
 * Then we are able to evaluate the deltaReLU and the second-layer multiplciation in 2 rounds.
 * The original sequential method is evaluate the deltaReLU firstly, costing 2 rounds, then evalute
 * the second layer multiplication, totally costing 3 rounds.
 * 
 * @param depth 
 * @param nCols 
 * @param origin 
 * @param maxIdx 
 * @return ShareBundle 
 */
ShareBundle ShareBundle::hierMaxpoolRowwise_opt(int depth, size_t nCols, const ShareBundle &origin, ShareBundle &maxIdx)const
{
    size_t ySize = nCols/2;
    if(nCols%2==1){ySize++;}
    int blkSize = 1<<depth;// depth starts from 1

    ShareBundle tmp(rows(), cols());
    tmp.shares = shares;
    ShareBundle blkBundle(rows(), nCols/2);// cond = 1, left > right -> max = (i*2)
    for(size_t i = 0; i < nCols/2; i++){
        blkBundle.shares.col(i) = shares.col(i*2) - shares.col(i*2+1);
    }

    ShareBundle cond(rows(), nCols/2);

    vector<ShareBundle>y;
    vector<BeaverTriple>triples;

    // This triple is for the max value.
    ShareBundle y0(rows(), nCols/2); //y0 = (a_{i*2} - a_{i*2+1})
    for(size_t i = 0; i < nCols/2; i++){
        // max(a_{i*2}, a_{i*2+1}) = cond_i * a_{i*2} + (1-cond_i) * a_{i*2+1}
        //                         = cond_i * (a_{i*2} - a_{i*2+1}) + a_{i*2+1}
        // So we can compute cond_i * (a_{i*2} - a_{i*2+1}) firstly.
        // But remember to add the degree-t term a_{i*2+1}
        y0.shares.col(i) = shares.col(i*2) - shares.col(i*2+1);
    }
    y.push_back(y0);
    triples.push_back(BeaverTriple(rows(), nCols/2));

    // These triples are for the max idx
    size_t leftBlkSize = blkSize/2;
    size_t rightBlkSize = (leftBlkSize*2 > origin.cols()) ? origin.cols() - leftBlkSize : leftBlkSize;
    if(rightBlkSize == leftBlkSize){// not the last layer
        for(size_t i = 0; i < leftBlkSize; i++){// store the leftBlk
            ShareBundle yi(rows(), nCols/2);
            for(size_t j = 0, k = i; j < nCols/2; j++, k+=blkSize){
                yi.shares.col(j) = maxIdx.shares.col(k);
            }
            y.push_back(yi);
            triples.push_back(BeaverTriple(rows(), nCols/2));
        }

        for(size_t i = 0; i < rightBlkSize; i++){// store the rightBlk
            ShareBundle yi(rows(), nCols/2);
            for(size_t j = 0, k = i+leftBlkSize; j < nCols/2; j++, k+=blkSize){
                yi.shares.col(j) = maxIdx.shares.col(k);
            }
            y.push_back(yi);
            triples.push_back(BeaverTriple(rows(), nCols/2));
        }
    }else{
        // In the last layer, for example size = 20. The leftBlkSize is 16 while the rightBlkSize is 4.
        // Then the size of cols of cond is only 1.
        assert(cond.shares.cols()==1);
        for(size_t i = 0; i < origin.cols(); i++){
            ShareBundle yi(rows(), 1);
            yi.shares = maxIdx.shares.col(i);

            y.push_back(yi);
            triples.push_back(BeaverTriple(rows(), 1));
        }
    }

    cond.shares = blkBundle.deltaReLU_opt(y, triples).shares;

    // Then decouple these triples to get the desired multiplication results.

    // First is the max value.
    ShareBundle maxBundle(rows(), ySize);
    maxBundle.shares.leftCols(nCols/2) = triples[0].mult().shares;
    for(size_t i = 0; i < nCols/2; i++){// Remember to add this value.
        maxBundle.shares.col(i) = maxBundle.shares.col(i) + shares.col(i*2+1);
    }
    if(nCols%2){
        maxBundle.shares.rightCols(1) = shares.rightCols(1);
    }
    
    size_t ttidx = 1;// triple idx
    if(rightBlkSize == leftBlkSize){// not the last layer
        for(size_t i = 0; i < leftBlkSize; i++){
            ShareBundle cond_idx = triples[ttidx].mult();//cond*y
            for(size_t j = 0, k = i; j < nCols/2; j++, k+=blkSize){
                maxIdx.shares.col(k) = cond_idx.shares.col(j);
            }
            ttidx++;
        }

        for(size_t i = 0; i < rightBlkSize; i++){// store the rightBlk
            ShareBundle cond_idx = triples[ttidx].x_times(-gfpScalar(1)).x_plus(1).mult();// (1-cond)*y
            for(size_t j = 0, k = i+leftBlkSize; j < nCols/2; j++, k+=blkSize){
                maxIdx.shares.col(k) = cond_idx.shares.col(j);
            }
            ttidx++;
        }
    }else{
        assert(cond.shares.cols()==1);
        for(size_t i = 0; i < leftBlkSize; i++){
            ShareBundle cond_idx = triples[ttidx].mult();//cond*y
            maxIdx.shares.col(i) = cond_idx.shares;
            ttidx++;
        }

        for(size_t j = leftBlkSize; j < origin.cols(); j++){
            ShareBundle cond_idx = triples[ttidx].x_times(-gfpScalar(1)).x_plus(1).mult();// (1-cond)*y
            maxIdx.shares.col(j) = cond_idx.shares;
            ttidx++;
        }
    }

    // Recursive termination point.
    if(ySize==1){
        return maxBundle;
    }

    return maxBundle.hierMaxpoolRowwise_opt(depth+1, ySize, origin, maxIdx);
}

// Only evaluate the maximum value on each row.
ShareBundle ShareBundle::hierMaxpoolRowwise(int depth, size_t nCols, const ShareBundle &origin)const
{
    size_t ySize = nCols/2;
    if(nCols%2==1){ySize++;}
    int blkSize = 1<<depth;// depth starts from 1

    ShareBundle blkBundle(rows(), nCols/2);
    for(size_t i = 0; i < nCols/2; i++){
        blkBundle.shares.col(i) = shares.col(i*2) - shares.col(i*2+1);
    }
    ShareBundle cond = blkBundle.deltaReLU(); // cond = 1, left > right -> max = (i*2)
    // cout<<cond.reveal()<<endl;

    ShareBundle maxBundle(rows(), ySize);
    // For verify
    gfpMatrix aTmp(rows(), nCols/2);
    gfpMatrix bTmp(rows(), nCols/2);
    gfpMatrix cTmp(rows(), nCols/2);
    for(size_t i = 0; i < nCols/2; i++){
        aTmp.col(i) = cond.shares.col(i);
        bTmp.col(i) = shares.col(i*2) - shares.col(i*2+1);
        maxBundle.shares.col(i) = cond.shares.col(i).array() * shares.col(i*2).array() + (1 - cond.shares.col(i).array()) * shares.col(i*2+1).array();
    }
    maxBundle.reduce_degree();
    // Verify
    for(size_t i = 0; i < nCols/2; i++){
        cTmp.col(i) = maxBundle.shares.col(i) - shares.col(i*2+1);
    }
    MultVerifier::push_triples(aTmp, bTmp, cTmp);
    aTmp.resize(0, 0);
    bTmp.resize(0, 0);
    cTmp.resize(0, 0);

    if(nCols%2==1){
        maxBundle.shares.rightCols(1) = shares.rightCols(1);
    }

    // Stop condition: Get the maximum.
    if(ySize==1){
        return maxBundle;
    }
    return maxBundle.hierMaxpoolRowwise(depth+1, ySize, origin);
}

// Only evaluates the maximum values.
ShareBundle ShareBundle::hierMaxpoolRowwise_opt(int depth, size_t nCols, const ShareBundle &origin)const
{
    size_t ySize = nCols/2;
    if(nCols%2==1){ySize++;}
    int blkSize = 1<<depth;// depth starts from 1


    ShareBundle blkBundle(rows(), nCols/2);// cond = 1, left > right -> max = (i*2)
    for(size_t i = 0; i < nCols/2; i++){
        blkBundle.shares.col(i) = shares.col(i*2) - shares.col(i*2+1);
    }

    ShareBundle cond(rows(), nCols/2);

    vector<ShareBundle>y;
    vector<BeaverTriple>triples;

    // This triple is for the max value.
    ShareBundle y0(rows(), nCols/2); //y0 = (a_{i*2} - a_{i*2+1})
    for(size_t i = 0; i < nCols/2; i++){
        // TODO: max(a_{i*2}, a_{i*2+1})=ReLU(a_{i*2+1}-a_{i*2}) + a_{i*2}
        // max(a_{i*2}, a_{i*2+1}) = cond_i * a_{i*2} + (1-cond_i) * a_{i*2+1}
        //                         = cond_i * (a_{i*2} - a_{i*2+1}) + a_{i*2+1}
        // So we can compute cond_i * (a_{i*2} - a_{i*2+1}) firstly.
        // But remember to add the degree-t term a_{i*2+1}
        y0.shares.col(i) = shares.col(i*2) - shares.col(i*2+1);
    }
    y.push_back(y0);
    triples.push_back(BeaverTriple(rows(), nCols/2));

    // Just delete the part for calculating the maxIdx

    cond.shares = blkBundle.deltaReLU_opt(y, triples).shares;

    // Then decouple these triples to get the desired multiplication results.

    // First is the max value.
    ShareBundle maxBundle(rows(), ySize);
    maxBundle.shares.leftCols(nCols/2) = triples[0].mult().shares;
    for(size_t i = 0; i < nCols/2; i++){// Remember to add this value.
        maxBundle.shares.col(i) = maxBundle.shares.col(i) + shares.col(i*2+1);
    }
    if(nCols%2){
        maxBundle.shares.rightCols(1) = shares.rightCols(1);
    }
    
    // Just delete the part for calculating the maxIdx

    // Recursive termination point.
    if(ySize==1){
        return maxBundle;
    }

    return maxBundle.hierMaxpoolRowwise_opt(depth+1, ySize, origin);
}

// Wrapper for Max Prime rowwise.
void ShareBundle::MaxPrimeRowwise(ShareBundle &maxIdx)const
{
    // *Sequential way: O(l)
    // ShareBundle maxRes(rows(), 1);
    // seqMaxpoolRowwise(maxRes, maxIdx);
    // return;

    // *Hierachical way: O(logl)
    // Remember to initalize.
    maxIdx.shares.setConstant(1);
    ShareBundle maxRes(rows(), 1);
    hierMaxpoolRowwise(1, cols(), *this, maxIdx);
    return;
}

// Wrapper for Max rowwise.
void ShareBundle::MaxRowwise(ShareBundle &maxRes)const
{
    maxRes = hierMaxpoolRowwise(1, cols(), *this);
    return;
}

// Wrapper for Max and Maxprime rowwise.
void ShareBundle::MaxRowwise(ShareBundle &maxRes, ShareBundle &maxIdx)const
{
    maxIdx.shares.setConstant(1);
    maxRes = hierMaxpoolRowwise(1, cols(), *this, maxIdx);
    return;
}

void ShareBundle::MaxRowwise_opt(ShareBundle &maxRes)const
{
    maxRes = hierMaxpoolRowwise_opt(1, cols(), *this);
    return;
}

void ShareBundle::MaxRowwise_opt(ShareBundle &maxRes, ShareBundle &maxIdx)const
{
    maxIdx.shares.setConstant(1);
    maxRes = hierMaxpoolRowwise_opt(1, cols(), *this, maxIdx);
    return;
}

// Bounding Power in FALCON.
// Bound power on each entry x to get alpha where 2^alpha <= x < 2^alpha+1 
gfpMatrix ShareBundle::bound_power()const
{
    // Not Used
}

gfpMatrix ShareBundle::bound_power_paralle() const
{
    // Not Used
}

gfpMatrix ShareBundle::bound_power_with_bits()const
{
    // Not Used
}
/********************************************************************************
 * 
 *       Definition of member functions about DoubleShareBundle
 * 
 * ******************************************************************************/

void DoubleShareBundle::input_from_random(int player_no)
{
    assert(degree==threshold);
    if(player_no == P->my_num()){
        random_matrix(secrets);
        octetStreams os(n_players);
        calculate_sharings(secrets, degree, shares,os);
        calculate_sharings(secrets, degree<<1, aux_shares, os);
        P->send_respective(os);
    }else{
        octetStream o;
        P->receive_player(player_no, o);
        unpack_matrix(shares, o);
        unpack_matrix(aux_shares, o);
    }
}

/**
 * @brief The request part the input_from_random in the multi-thread version.
 * The distributer calculates the sharings of degree t and of degree 2t, and only requests to send them.
 * The receiver only requests to receive the shares.
 * 
 * @param player_no 
 * @param os_send 
 * @param o_receive 
 */
void DoubleShareBundle::input_from_random_request(int player_no, octetStreams & os_send, octetStream & o_receive)
{
    if(player_no == P->my_num()){
        random_matrix(secrets);
        os_send.reset(n_players);
        calculate_sharings(secrets, degree, shares, os_send);
        calculate_sharings(secrets, degree<<1, aux_shares, os_send);
        P->request_send_respective(os_send);
    }else{P->request_receive(player_no, o_receive);}
}

/**
 * @brief The wait part of the multi-thread version for inputting sharings.
 * The distribtuer waits to send the shares(and aux_shares) to the corresponding party.
 * The receiver waits to receive the shares(and aux_shares) from the distributer.
 * @param player_no 
 * @param os_send 
 * @param o_receive 
 */
void DoubleShareBundle::finish_input_from(int player_no, octetStreams &os_send, octetStream &o_receive)
{
    if(player_no == P->my_num()){
        P->wait_send_respective(os_send);
    }else{
        P->wait_receive(player_no, o_receive);
        unpack_matrix(shares, o_receive);
        unpack_matrix(aux_shares, o_receive);
    }
}

// * With PRG
void DoubleShareBundle::input_from_random_PRG(int player_no)
{
    assert(degree==threshold);
    if(player_no==P->my_num()){
        random_matrix(secrets);
        octetStreams os(n_party_PRG());
        calculate_t_sharings_PRG(secrets, shares, os);
        P->send_respective(start_party_PRG(), n_party_PRG(), os);
        calculate_2t_sharings_PRG(secrets, aux_shares);
    }else{
        get_t_sharings_PRG(player_no, shares);
        get_2t_sharings_PRG(player_no, aux_shares);
    }
}

// 1. Random secrets and calculate the t-sharings, and request to send.
// void DoubleShareBundle::input_from_random_request_PRG(int player_no, octetStreams &os_send, gfpMatrix &shares_prng, octetStream &o_receive)
// {
//     assert(degree==threshold);
//     if(player_no == P->my_num()){
//         random_matrix(secrets);
//         os_send.reset(n_party_PRG());
//         shares_prng.resize(degree, size());

//         calculate_t_sharings_PRG(secrets, shares, os_send);
//         P->request_send_respective(start_party_PRG(), n_party_PRG(), os_send);
//     }else{
//         shares_prng.resize(degree, size());
//         get_t_sharings_PRG_request(player_no, shares_prng, o_receive);
//     }
// }

// 2. Calculate the 2t-sharings without communication. 
void DoubleShareBundle::input_aux_from_random_PRG(int player_no)
{
    if(player_no == P->my_num()){
        calculate_2t_sharings_PRG(secrets, aux_shares);
    }else{
        get_2t_sharings_PRG(player_no, aux_shares);
    }
}

// 3. Wait to receive the t-sharing.
// void DoubleShareBundle::finish_input_from_PRG(int player_no, octetStreams &os_send, const gfpMatrix &shares_prng, octetStream &o_receive)
// {
//     if(player_no == P->my_num()){
//         P->wait_send_respective(start_party_PRG(), n_party_PRG(), os_send);
//     }else{
//         get_t_sharings_PRG_wait(player_no, shares_prng, shares, o_receive);
//     }
// }

gfpMatrix DoubleShareBundle::reveal_aux(int degree)
{
    ShareBundle aux(aux_shares.rows(), aux_shares.cols());
    aux.shares = aux_shares;
    aux.set_degree(degree);
    return aux.reveal(true);
}

DoubleShareBundle& DoubleShareBundle::reduced_random()
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline()){
            Phase->generate_reduced_random_sharings(size());
        }
        else{
            Phase->switch_to_offline();
            Phase->generate_reduced_random_sharings(size());
            Phase->switch_to_online();
        }
    }
    // ( [r]_t, [r]_2t )
    DoubleRandom::get_random_pairs(DoubleRandom::queueReducedRandom, shares, aux_shares);
    return *this;
}

DoubleShareBundle& DoubleShareBundle::truncated_random(ShareBundle &msb)
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    msb.shares.setConstant(0);
    return *this;
#endif
    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_truncated_random_sharings(size());
        else{
            Phase->switch_to_offline();
            Phase->generate_truncated_random_sharings(size());
            Phase->switch_to_online();
        }
    }
    // ( [r/2^d]_t, [r]_t )
    // DoubleRandom::get_random_pairs(DoubleRandom::queueTruncatedRandom, shares, aux_shares);
    DoubleRandom::get_random_triples(DoubleRandom::queueTruncatedRandom, shares, aux_shares, msb.shares);
    return *this;
}

DoubleShareBundle& DoubleShareBundle::truncated_random(size_t precision)
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_truncated_random_sharings(size(), precision);
        else{
            Phase->switch_to_offline();
            Phase->generate_truncated_random_sharings(size(), precision);
            Phase->switch_to_online();
        }
    }
    // ( [r/2^p]_t, [r]_t )
    DoubleRandom::get_random_pairs(DoubleRandom::queueTruncatedRandomInML, shares, aux_shares);
    return *this;
}

DoubleShareBundle& DoubleShareBundle::reduced_truncated_random(ShareBundle &msb)
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    aux_aux_shares.setConstant(0);
    msb.shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_reduced_truncated_sharings(size());
        else{
            Phase->switch_to_offline();
            Phase->generate_reduced_truncated_sharings(size());
            Phase->switch_to_online();
        }
    }
    // ( [r/2^d]_t, [r]_2t, [r]_t [r_msb] )
    DoubleRandom::get_random_quadruples(DoubleRandom::queueReducedTruncatedRandom, shares, aux_shares, aux_aux_shares, msb.shares);
    return *this;
}

DoubleShareBundle& DoubleShareBundle::reduced_truncated_random(size_t logLearningRate, size_t logMiniBatch)
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    aux_aux_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_reduced_truncated_sharings(size(), logLearningRate, logMiniBatch);
        else{
            Phase->switch_to_offline();
            Phase->generate_reduced_truncated_sharings(size(), logLearningRate, logMiniBatch);
            Phase->switch_to_online();
        }
    }
    // ( [r/2^d]_2t, [r]_2t, [r]_t)
    DoubleRandom::get_random_triples(DoubleRandom::queueReducedTruncatedInML, shares, aux_shares, aux_aux_shares);
    return *this;
}

DoubleShareBundle& DoubleShareBundle::reduced_truncated_random(size_t precision)
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    aux_aux_shares.setConstant(0);
    return *this;
#endif

    // * We explicitly generate the reduced and truncated random sharings in the functionality (switch to offline).
    // if(!Phase->is_true_offline()){
    //     if(Phase->is_Offline())
    //         Phase->generate_reduced_truncated_sharings(size(), precision);
    //     else{
    //         Phase->switch_to_offline();
    //         Phase->generate_reduced_truncated_sharings(size(), precision);
    //         Phase->switch_to_online();
    //     }
    // }
    // ( [r/2^precision]_2t, [r]_2t, [r]_t)
    DoubleRandom::get_random_triples(DoubleRandom::queueReducedTruncatedWithPrecisionRandom, shares, aux_shares, aux_aux_shares);
    return *this;
}

DoubleShareBundle& DoubleShareBundle::reduced_truncated_random(vector<size_t> &precision)
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    aux_aux_shares.setConstant(0);
    return *this;
#endif

    // * We explicitly generate the reduced and truncated random sharings in the functionality (switch to offline).
    DoubleRandom::get_random_triples(DoubleRandom::queueReducedTruncatedWithPrecisionRandom, shares, aux_shares, aux_aux_shares);
    return *this;
}

DoubleShareBundle& DoubleShareBundle::unbounded_prefix_mult_random()
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_unbounded_mult_random_sharings(rows(), cols());
        else{
            Phase->switch_to_offline();
            Phase->generate_unbounded_mult_random_sharings(rows(), cols());
            Phase->switch_to_online();
        }
    }
    // Each row is an instance of unbounded prefix mult.
    //  ( [b1]_t , [b1^-1]_t )
    //  ( [b2]_t , [b1 * b2^-1]_t,)
    //  ( [b3]_t , [b2 * b3^-1]_t)
    //  ...
    //  ( [bi]_t , [bi-1 * bi^-1]_t) for i = 2, ..., l
    DoubleRandom::get_random_pairs(DoubleRandom::queueUnboundedMultRandom, shares, aux_shares);
    return *this;
}

DoubleShareBundlePair& DoubleShareBundlePair::mult_reduced_random_pair()
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    other_shares.setConstant(0);
    other_aux_shares.setConstant(0);
    prod_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_mult_reduced_random_sharing_pairs(size());
        else{
            Phase->switch_to_offline();
            Phase->generate_mult_reduced_random_sharing_pairs(size());
            Phase->switch_to_online();
        }
    }

    DoubleRandom::get_mult_reduced_random_sharing_pairs(shares, aux_shares, other_shares, other_aux_shares, prod_shares);
    return *this;
}

ReusedDoubleShareBundlePair& ReusedDoubleShareBundlePair::reused_mult_reduced_random_pairs()
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    other_shares.setConstant(0);
    other_aux_shares.setConstant(0);
    prod_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_reused_mult_reduced_random_sharing_pairs(size());
        else{
            Phase->switch_to_offline();
            Phase->generate_reused_mult_reduced_random_sharing_pairs(size());
            Phase->switch_to_online();
        }
    }

    DoubleRandom::get_reused_mult_reduced_random_sharing_pairs(shares, aux_shares, other_shares, other_aux_shares, prod_shares);
    return *this;
}

ReusedDoubleShareBundlePair16& ReusedDoubleShareBundlePair16::reused_mult_reduced_random_pairs_16()
{
#ifdef ZERO_OFFLINE
    shares.setConstant(0);
    aux_shares.setConstant(0);
    other_shares.setConstant(0);
    other_aux_shares.setConstant(0);
    prod_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_reused_mult_reduced_random_sharing_pairs_16(size());
        else{
            Phase->switch_to_offline();
            Phase->generate_reused_mult_reduced_random_sharing_pairs_16(size());
            Phase->switch_to_online();
        }
    }

    DoubleRandom::get_reused_mult_reduced_random_sharing_pairs_16(shares, aux_shares, other_shares, other_aux_shares, prod_shares);
    return *this;
}

BeaverBundleTriple& BeaverBundleTriple::beaver_triple(bool is_offline)
{
#ifdef ZERO_OFFLINE
    left_shares.setConstant(0);
    right_shares.setConstant(0);
    prod_shares.setConstant(0);
    return *this;
#endif

    if(!Phase->is_true_offline()){
        if(Phase->is_Offline())
            Phase->generate_beaver_triple(x_size(), y_size(), z_size(), is_offline);
        else{
            Phase->switch_to_offline();
            Phase->generate_beaver_triple(x_size(), y_size(), z_size(), is_offline);
            Phase->switch_to_online();
        }
    }

    DoubleRandom::get_beaver_triple(left_shares, right_shares, prod_shares);
    return *this;
}

gfpMatrix BeaverBundleTriple::beaver_mult(const gfpMatrix &input_x, const gfpMatrix &input_y)
{
    assert(x_size() == input_x.rows());
    assert(y_size() == input_x.cols());
    assert(y_size() == input_y.rows());
    assert(z_size() == input_y.cols());

    ShareBundle zShareBundle(x_size()*y_size()+y_size()*z_size(), 1);
    zShareBundle.shares(seqN(0, x_size()*y_size()), 0) = (input_x-left_shares).reshaped<RowMajor>(); // u = [x]_t - [a]_t
    zShareBundle.shares(seqN(x_size()*y_size(), y_size()*z_size()), 0) = (input_y - right_shares).reshaped<RowMajor>(); // v = [y]_t - [b]_t
    zShareBundle.reveal(true);

    gfpMatrix u = zShareBundle.secret()(seqN(0, x_size()*y_size()), 0).reshaped<RowMajor>(x_size(), y_size());
    gfpMatrix v = zShareBundle.secret()(seqN(x_size()*y_size(), y_size()*z_size()), 0).reshaped<RowMajor>(y_size(), z_size());
    gfpMatrix res = prod_shares 
                    + u * right_shares
                    + left_shares * v
                    + u * v;
    return res;
}


}
