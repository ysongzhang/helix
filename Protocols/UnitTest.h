#ifndef PROTOCOLS_UNITTEST_H_
#define PROTOCOLS_UNITTEST_H_
#include "Protocols/PhaseConfig.h"
using namespace std;
namespace hmmpc
{
// Share with PRG
void debugSharePRG(PhaseConfig *phase);
void debugShareBundlePRG(PhaseConfig *phase);
void debugRandomSharePRG(PhaseConfig *phase);

// ShareBundle
void debugUnboundedPrefixMult(PhaseConfig *phase);
void debugGetLSB(PhaseConfig *phase);
void debugGetMSB(PhaseConfig *phase);
void debugReLU(PhaseConfig *phase);
void debugBoundPower(PhaseConfig *phase);
void debugMaxPool(PhaseConfig *phase);
void debug3LayerMult(PhaseConfig *phase);
void debug3LayerMultReused(PhaseConfig *phase);
void debugUnboundedMult(PhaseConfig *phase);
void debugPrefixMultDN(PhaseConfig *phase);

// Bit
void debugBitOp(PhaseConfig *phase);
void debugUnboundedOp(PhaseConfig *phase);
void debugPrefixBlkOp(PhaseConfig *phase);
void debugPrefixOp(PhaseConfig *phase);

void debugLessThanUnsigned(PhaseConfig *phase);
void debugBitAdd(PhaseConfig *phase);
void debugBitsDecomposition(PhaseConfig *phase);

void debugCarrayPropagation(PhaseConfig *phase);

// Bench
void benchReLU(PhaseConfig *phase, int num, int true_offline);

// Debug
void debugGeneral(PhaseConfig *phase);
}

#endif