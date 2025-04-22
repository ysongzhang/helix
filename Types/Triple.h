#ifndef TYPES_TRIPLE_H_
#define TYPES_TRIPLE_H_

#include "Math/gfpScalar.h"
#include "Math/gfpMatrix.h"

namespace hmmpc
{

// The product could be t-sharing or just public value. Both cases have the same procedure.
class Triple;  // gfpScalar class
class DotTriple;
class MatrixTriple;  // gfpMatrix class
class MatrixDimension;


class Triple 
{
public:
    gfpScalar left;
    gfpScalar right;
    gfpScalar product;  // product = left * right 

    Triple(const gfpScalar a, const gfpScalar b, const gfpScalar c): left(a), right(b), product(c){}    
    Triple(const gfpScalar a, const gfpScalar b): left(a), right(b){product = 0;}
};

// TODO: may be using reference
class DotTriple
{
public:
    gfpVector left;  // gfpVector
    gfpVector right;
    gfpScalar product;  // product = left * right

    // Construction
    DotTriple(){}
    DotTriple(const gfpVector a, const gfpVector b, const gfpScalar c): left(a), right(b), product(c){}
};


// Matrix multiplication
// TODO: may be using reference
class MatrixTriple
{
public:
    gfpMatrix left;
    gfpMatrix right;
    gfpMatrix product;  // product = left * right

    // Construction
    MatrixTriple(){}
    MatrixTriple(const gfpMatrix a, const gfpMatrix b): left(a), right(b){}
    MatrixTriple(const gfpMatrix a, const gfpMatrix b, const gfpMatrix c): left(a), right(b), product(c){
        assert(a.rows() == c.rows());
        assert(a.cols() == b.rows());
        assert(b.cols() == c.cols());
    }

    // Phrase to Dot product
    vector<DotTriple> phrase_to_dot_products();
    
    // Get dimensions 
    size_t first_dimension(){return left.rows();}
    size_t second_dimension(){return left.cols();}
    size_t third_dimension(){return right.cols();}

    // Compress algorithm with random elements
    void compress(const gfpVector &randomness_left, const gfpVector &randomness_right, gfpVector &compressed_left, gfpVector &compressed_right, gfpScalar &compressed_prod);
};

class MatrixDimension
{
public:
    size_t dim_x;
    size_t dim_y;
    size_t dim_z;

    // Construction
    MatrixDimension(){}
    MatrixDimension(const size_t a, const size_t b, const size_t c): dim_x(a), dim_y(b), dim_z(c){}

    size_t both_dimension(){return dim_x*dim_y+dim_y*dim_z;}
    size_t z_dimension(){return dim_x*dim_z;}
};

// LOG: note that the class function defined in .h files must be inline.
inline vector<DotTriple> MatrixTriple::phrase_to_dot_products()
{
    vector<DotTriple> res;
    DotTriple tb;
    gfpVector a(left.cols());
    gfpVector b(right.rows());
    gfpScalar c;
    for(size_t i = 0; i < product.rows(); i++) {
        for (size_t j = 0; j < product.cols(); j++)
        {
            a = left.row(i).transpose();
            b = right.col(j);
            c = product(i, j);
            tb = DotTriple(a, b, c);
            res.push_back(tb);
        }
    }
    return res;
}

inline void MatrixTriple::compress(const gfpVector &randomness_left, const gfpVector &randomness_right, gfpVector &compressed_left, gfpVector &compressed_right, gfpScalar &compressed_prod)
{
    assert(first_dimension() == randomness_left.size());
    assert(third_dimension() == randomness_right.size());
    assert(second_dimension() == compressed_left.size() && second_dimension() == compressed_right.size());
    compressed_left.setConstant(0);
    compressed_right.setConstant(0);
    for (size_t i = 0; i < first_dimension(); i++)
    {
        compressed_left += (randomness_left(i) * left.row(i)).transpose();
    }
    for (size_t i = 0; i < third_dimension(); i++)
    {
        compressed_right += (randomness_right(i) * right.col(i));
    }
    gfpMatrix randomness_prod = randomness_left.rowwise().replicate(third_dimension());
    randomness_prod = randomness_prod.array() * randomness_right.transpose().colwise().replicate(first_dimension()).array();
    compressed_prod = (randomness_prod.array() * product.array()).sum();
}

}
#endif