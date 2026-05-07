#include "Math/UnitTest.h"

using namespace std;
namespace hmmpc
{
void testGfp()
{
    cout<<"[UnitTest for Gfp]:"<<endl;

    gfpScalar x=2;
    cout<<"x = "<<x<<endl;
    gfpScalar minus_x = -x;
    gfpScalar inverse_x = multiplicative_inverse((TYPE)(2));

    cout<<"-x = "<<minus_x<<endl;
    cout<<"x^{-1} = "<<inverse_x<<endl;
    cout<<"x+(-x)=" <<x+minus_x<<endl;
    cout<<"x*(x^-1) = "<<x*inverse_x<<endl;
}

void debugGfpDivision()
{
    // Not Used
}

void debugCNNExtend()
{
    size_t B = 1;
    size_t iw = 4;
    size_t ih = 4;
    size_t Din = 1;
    size_t P = 2;
    size_t f = 2;
    size_t S = 2;
    size_t Dout = 1;
    size_t ow 	= (((iw-f+2*P)/S)+1);
	size_t oh	= (((ih-f+2*P)/S)+1);
    gfpMatrix input(B,iw*ih*Din);
    gfpMatrix paddedInput(B, (iw+2*P)*(ih+2*P)*Din);
    gfpMatrix extendInput(B*oh*ow, f*f*Din);

    input<<1,2,3,4,5,6,7,8,9,10,
    11,12,13,14,15,16;
    cout<<input.reshaped<RowMajor>(iw, ih)<<endl;
    zeroPad(input, paddedInput, iw, ih, P, Din, B);
    cout<<paddedInput.reshaped<RowMajor>(iw+2*P, ih+2*P)<<endl;

    convolExtend(paddedInput, extendInput, iw, ih, ow, oh, Din, S, f, P, B);
    cout<<extendInput<<endl;
}
void debugMaxpoolExtend()
{
    gfpMatrix A(1,9*2);
    A<<1,2,3,4,5,6,7,8,9,
    11,12,13,14,15,16,17,18,19;
    cout<<A.reshaped<RowMajor>(6, 3)<<endl;
    
    gfpMatrix C(8, 4);
    maxpoolExtend(A, C, 3, 3, 2, 2, 2, 1, 2, 1);
    cout<<C<<endl;
}
}