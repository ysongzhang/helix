#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include "Types/sfix.h"
#include "Types/sintMatrix.h"
#include "Types/cfixMatrix.h"
#include "Types/sfixMatrix.h"
#include "Protocols/UnitTest.h"
#include "Types/UnitTest.h"

using namespace std;
using namespace hmmpc;
int main(int argc, char** argv)
{

    int player_no = atoi(argv[1]);
    int portnum_base = 6000;
    string filename = "Test/HOSTS.example";
    Names player_name = Names(player_no, portnum_base, filename);
    ThreadPlayer P(player_name, "");
    int threshold = 1;

    PhaseConfig phase;
    phase.init(P.num_players(), threshold, &P);

    string which_test = argv[2];
    if (which_test.compare("Trunc")==0)
    {
        // Prot 3.1: Multiplication with secret fixed-point number and public fixed-point number
        debugSfixMatMulCfix(&phase);
    }
    else if (which_test.compare("Fixed-Mult")==0)
    {
        // Prot 3.2: Multiplication with secret fixed-point number and secret fixed-point number
        debugSfixMatMul(&phase);
    }
    else if (which_test.compare("PreOR")==0)
    {
        // Prot 4.1: Compute prefix-OR of secret shared bits
        debugPrefixOp(&phase);
    }
    else if (which_test.compare("Bitwise-LT")==0)
    {
        // Prot 4.2: Bitwise-LessThan protocol
        debugLessThanUnsigned(&phase);
    }
    else if (which_test.compare("DReLU")==0)
    {
        // Prot 5.1: DReLU
        debugGetMSB(&phase);
    }
    else if (which_test.compare("ReLU")==0)
    {
        // Prot 5.3: ReLU
        debugReLU(&phase);
    }
    else if (which_test.compare("Maxpool")==0)
    {
        // Prot 5.4: ReLU
        debugMaxPool(&phase);
    }
    else if (which_test.compare("3L-DN")==0)
    {
        debug3LayerMult(&phase);
    }
    else if (which_test.compare("3L-DN-reused")==0)
    {
        debug3LayerMultReused(&phase);
    }
    else if (which_test.compare("Unbounded-mult")==0)
    {
        debugUnboundedMult(&phase);
    }
    else if (which_test.compare("Prefix-mult")==0)
    {
        debugPrefixMultDN(&phase);
    }
    else if (which_test.compare("Bench-ReLU")==0)
    {
        benchReLU(&phase, atoi(argv[3]), atoi(argv[4]));
    }
    else
    {
        // *Other unit tests;

        // debugSintMatMul(&phase);
        // debugUnboundedPrefixMult(&phase);

        // testSfixMul(&phase);
        // testSintMul(&phase);
        // testCint(&phase);
        // testSfix(&phase);

        // debugBitOp(&phase);
        // debugGetLSB(&phase);        

        // debugSharePRG(&phase);
        // debugShareBundlePRG(&phase);
        // debugRandomSharePRG(&phase);
        debugGeneral(&phase);
    }
    phase.print_offline_communication();
    phase.print_online_communication();

    return 0;
}
