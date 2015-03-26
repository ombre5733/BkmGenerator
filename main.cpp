// The MIT License (MIT)
//
// Copyright (c) 2015 Manuel Freiberger
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>

using namespace std;


constexpr long double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899l;

unsigned clog2(unsigned long long value)
{
    unsigned long long log = 0;
    while ((1ull << log) < value)
        ++log;
    return log;
}

struct Config
{
    Config(int numBits, int numFractionBits, long long ut)
        : NB(numBits),
          FB(numFractionBits),
          underflowThreshold(ut)
    {
        if (FB < 4)
            throw runtime_error("Need at least 4 fraction bits.");
        if (NB < FB + 1)
            throw runtime_error("Must have more total bits than fraction bits.");
    }

    string dataType() const
    {
        stringstream ss;
        ss << "bit signed [" << HB << ":" << -FB << "]";
        return ss.str();
    }

    string format(long long v) const
    {
        stringstream ss;
        if (v < 0)
            ss << '-';
        else
            ss << ' ';
        ss << NB << "'sd";
        ss << (v < 0 ? -v : v);
        return ss.str();
    }

    long long fixed(long double v) const
    {
        long long scal = 1ll << FB;
        return static_cast<long long>(std::floor(v * scal + 0.5)); // TODO: correct the rounding
    }

    long long ceilToFixed(long double v) const
    {
        if (v < 0)
            throw runtime_error("Must be positive");

        long long scal = 1ll << FB;
        return static_cast<long long>(std::ceil(v * scal));
    }

    long long floorToFixed(long double v) const
    {
        if (v < 0)
            throw runtime_error("Must be positive");

        long long scal = 1ll << FB;
        return static_cast<long long>(std::floor(v * scal));
    }

    long long maxFixed() const
    {
        return (1ll << (NB - 1)) - 1;
    }


    // Total number of bits.
    int NB;
    // Number of fraction bits
    int FB;

    // Number of integer bits.
    const int IB = NB - FB;
    // Highest bit index.
    const int HB = NB - FB - 1;

    // The maximum number of shifts which will be applied in the reduction phase.
    int maxNumReductionShifts;

    // If the magnitude of the input is smaller than this value,
    // no logarithm is computed.
    const long long underflowThreshold;
};

struct ConfigBuilder
{
    ConfigBuilder()
        : totalNumBits(16),
          numFractionBits(12),
          underflowThreshold(1)
    {
    }

    ConfigBuilder& setNumBits(unsigned b)
    {
        totalNumBits = b;
        return *this;
    }

    ConfigBuilder& setNumFractionBits(unsigned b)
    {
        numFractionBits = b;
        return *this;
    }

    ConfigBuilder& setUnderflowThreshold(long long t)
    {
        underflowThreshold = t;
        return *this;
    }

    Config create()
    {
        return Config(totalNumBits, numFractionBits, underflowThreshold);
    }

    unsigned totalNumBits;
    unsigned numFractionBits;
    long long underflowThreshold;
};

/*
 * Range reduction:
 *
 * if (x < 0) x := -x,    y := -y
 * if (y < 0) x := -y,    y := x
 * if (y > x) x := x - y, y := x + y
 *
 * (a + ib) * e^(i pi)    = (a + ib) * -1 = -a - ib
 * (a + ib) * e^(i pi/2)  = (a + ib) * i  = -b + ia
 * (a + ib) * e^(-i pi/4) = (a + ib) * (1 + i) / sqrt(2) = (a - b + ia + ib) / sqrt(2)
*/

class Counter
{
public:
    Counter(string _name, unsigned _numBits)
        : name(_name),
          numBits(_numBits)
    {
    }

    static Counter fromSize(string name, unsigned size)
    {
        unsigned numBits = clog2(size);
        stringstream ss;
        ss << "r" << numBits << name;
        return Counter(ss.str(), numBits);
    }

    static Counter fromMaxValue(string name, unsigned long long maxValue)
    {
        unsigned numBits = clog2(maxValue + 1);
        stringstream ss;
        ss << "r" << numBits << name;
        return Counter(ss.str(), numBits);
    }

    string name;
    unsigned numBits;
};

void lModeComputer(std::ofstream& file, const Config& cfg)
{
    const int numIter = cfg.FB;
    const int bkmLookupSize = 8 * numIter;

    // The convergence bounds.
    long long lowerConvergenceBound = cfg.ceilToFixed(0.64);
    long long upperConvergenceBound = cfg.floorToFixed(1.4);
    if (lowerConvergenceBound > cfg.maxFixed())
        throw runtime_error("Cannot represent the value 0.64 in the chosen fixed point format.");
    if (upperConvergenceBound <= 0)
        throw runtime_error("Cannot represent the value 1.4 in the chosen fixed point format.");

    // Determine how often we have to shift the largest and smallest positive
    // value until it is within the convergence bounds.
    auto isInBounds = [&](long long v) -> bool {
        return v >= lowerConvergenceBound && v <= upperConvergenceBound;
    };
    int requiredLeftShifts = 0;
    {
        long long min = 1;
        while (!isInBounds(min))
        {
            min <<= 1;
            ++requiredLeftShifts;
        }
    }
    int requiredRightShifts = 0;
    {
        long long max = cfg.maxFixed();
        while (!isInBounds(max))
        {
            max >>= 1;
            ++requiredRightShifts;
        }
    }
    int requiredShifts = max(requiredLeftShifts, requiredRightShifts);

    long long underflowThreshold = lowerConvergenceBound >> requiredShifts;
    if (underflowThreshold == 0)
        underflowThreshold = 1;
    long long overflowThreshold = upperConvergenceBound << requiredShifts;
    if (overflowThreshold >= cfg.maxFixed())
        overflowThreshold = -1;




    Counter iterationCounter = Counter::fromMaxValue("IterationCounter", numIter);
    Counter lookupOffset = Counter::fromSize("LookupOffset", bkmLookupSize);
    Counter reductionShiftCounter = Counter::fromMaxValue("ReductionShifts", requiredShifts);

#if 0
    if (cfg.underflowThreshold > upperConvergenceBound)
        throw runtime_error("The smallest magnitude is too large.");

    // Determine x such that the smallest possible magnitude times 2^x is
    // in the convergence bound.
    int maxScalingLeftShifts = 0;
    {
        long long t = cfg.underflowThreshold;
        while (t < lowerConvergenceBound)
        {
            t <<= 1;
            ++maxScalingLeftShifts;
        }
    }

    // In the worst case, we have to shift by as many places as there are
    // fraction bits or integer bits. If all bits are set, we have to
    // shift until the integer part is zero. If only the LSB is set, we have
    // to shift until the fractional part is zero.
    // More mathematically, we scale by 2^x, with -IB <= x <= FB.
    // We encode the scaling as 2^(y-b), with b = IB and 0 <= y <= IB + FB = NB.
    const int reductionShiftCounter.numBits = clog2(cfg.NB + 1);
    string reductionShiftCounter.name;
    {
        stringstream ss;
        ss << "r" << reductionShiftCounter.numBits << "Scaling";
        reductionShiftCounter.name = ss.str();
    }
#endif


    file << "// The MIT License (MIT)\n"
            "//\n"
            "// Copyright (c) 2015 Manuel Freiberger\n"
            "//\n"
            "// Permission is hereby granted, free of charge, to any person obtaining a copy\n"
            "// of this software and associated documentation files (the \"Software\"), to deal\n"
            "// in the Software without restriction, including without limitation the rights\n"
            "// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
            "// copies of the Software, and to permit persons to whom the Software is\n"
            "// furnished to do so, subject to the following conditions:\n"
            "//\n"
            "// The above copyright notice and this permission notice shall be included in all\n"
            "// copies or substantial portions of the Software.\n"
            "//\n"
            "// THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
            "// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
            "// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
            "// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
            "// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
            "// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
            "// SOFTWARE.\n\n\n";

    file << "`timescale 1ns/1ns\n\n";
    file << "`ifndef AssignDelay\n"
         << "    `define AssignDelay #1\n"
         << "`endif\n\n";

    file << "module BkmLComputer\n"
         << "       (input bit clock,\n"
         << "        input bit reset_n,\n"
         << "        input " << cfg.dataType() << " L_in_r,\n"
         << "        input " << cfg.dataType() << " L_in_i,\n"
         << "        input " << cfg.dataType() << " E_in_r,\n"
         << "        input " << cfg.dataType() << " E_in_i,\n"
         << "        input bit inputValid,\n"
         << "        output " << cfg.dataType() << " L_out_r,\n"
         << "        output " << cfg.dataType() << " L_out_i,\n"
         << "        output " << cfg.dataType() << " E_out_r,\n"
         << "        output " << cfg.dataType() << " E_out_i,\n"
         << "        output bit outputValid,\n"
         << "        output bit outOfRange\n"
         << "       );\n\n";

    file << "    // L-mode algorithm\n"
         << "    //     E_1 = E\n"
         << "    //     L_1 = 2 (L - 1)\n"
         << "    //     E_{k+1} = E_k - ln(1 + d_k 2^{-k})\n"
         << "    //     L_{k+1} = 2 (L_k + d_k) + 2^{-k+1} L_k d_k\n"
         << "    //\n"
         << "    //     Stage 1:\n"
         << "    //         d_k = f(L_k)\n"
         << "    //     Stage 2:\n"
         << "    //         Delta E_k^r = LUT^r[k, d_k]\n"
         << "    //         H_k = L_k * d_k\n"
         << "    //         L_k = L_k + d_k\n"
         << "    //     Stage 3:\n"
         << "    //         Delta E_k^i = LUT^i[k, d_k]\n"
         << "    //         H_k = H_k >>> (k-1)\n"
         << "    //         L_k = 2 * L_k\n"
         << "    //     Stage 4:\n"
         << "    //         E_k = E_k - Delta E_k\n"
         << "    //         L_k = L_k + H_k\n"
         << "    //\n"
         << "    //     Result:\n"
         << "    //         E_k -> E + ln(L)\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     Constants\n";
    file << "    // ----============================================================----\n\n";

    file << "    // Convergence bounds of the real part.\n"
         << "    localparam " << cfg.dataType() << " CONST_0_64 = " << cfg.format(lowerConvergenceBound) << ";\n"
         << "    localparam " << cfg.dataType() << " CONST_1_40 = " << cfg.format(upperConvergenceBound) << ";\n\n";

    file << "    // The bounds which we have to process during the scaling stage in the reduction phase.\n"
         << "    localparam " << cfg.dataType() << " CONST_UNDERFLOW_THRESHOLD = "
         << cfg.format(underflowThreshold) << ";\n";
    if (overflowThreshold > 0)
    {
         file << "    localparam " << cfg.dataType() << " CONST_OVERFLOW_THRESHOLD = "
              << cfg.format(overflowThreshold) << ";\n";
    }
    file << "\n";

    file << "    // Pi and fractions thereof.\n"
         << "    localparam " << cfg.dataType() << " CONST_PI   = " << cfg.format(cfg.fixed(pi)) << ";\n"
         << "    localparam " << cfg.dataType() << " CONST_PI_2 = " << cfg.format(cfg.fixed(pi / 2)) << ";\n"
         << "    localparam " << cfg.dataType() << " CONST_PI_4 = " << cfg.format(cfg.fixed(pi / 4)) << ";\n\n";

    file << "    // log(1/sqrt(2))\n"
         << "    localparam " << cfg.dataType() << " CONST_LOG_1_SQRT2 = "
         << cfg.format(cfg.fixed(log(1.0 / sqrt(2.0)))) << ";\n"
         << "    // log(2/sqrt(5))\n"
         << "    localparam " << cfg.dataType() << " CONST_LOG_2_SQRT5 = "
         << cfg.format(cfg.fixed(log(2.0 / sqrt(5.0)))) << ";\n"
         << "    // atan(1/2)\n"
         << "    localparam "<< cfg.dataType() << " CONST_ATAN_1_2 = "
         << cfg.format(cfg.fixed(atan(1.0 / 2.0))) << ";\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     Lookup table\n";
    file << "    // ----============================================================----\n\n";

    auto rightJustify = [&] (string s) -> string {
        stringstream ss;
        ss << setw(std::ceil(log10(cfg.maxFixed())) + 6) << s;
        return ss.str();
    };
    auto log = [&] (unsigned i, long double dr, long double di) -> string {
        long long pow = 1ll << i;
        auto x = 1.0 + dr / pow;
        auto y = di / pow;
        return rightJustify(cfg.format(cfg.fixed(0.5 * std::log(x*x + y*y))));
    };

    auto arctan = [&] (unsigned i, long double dr) -> string {
        long long pow = 1ll << i;
        return rightJustify(cfg.format(cfg.fixed(std::atan(1.0 / (pow + dr)))));
    };

    file << "    " << cfg.dataType() << " LUT_BKM[" << bkmLookupSize << "] = '{";
    for (int i = 1; i <= numIter; ++i)
    {
        if (i > 1)
            file << ",";
        file << "\n        ";
        file << log(i, -1, 0) << ", ";
        file << log(i, -1, 1) << ", ";
        file << log(i,  0, 1) << ", ";
        file << log(i,  1, 0) << ", ";
        file << log(i,  1, 1) << ", ";
        file << arctan(i, -1) << ", ";
        file << arctan(i,  0) << ", ";
        file << arctan(i,  1);
    }
    file << "};\n\n";

    file << "    " << cfg.dataType() << " LUT_LOG_POW2[" << requiredShifts + 1 << "] = '{";
    for (int i = 0; i <= requiredShifts; ++i)
    {
        if (i != 0)
            file << ",";
        file << "\n        ";
        file << rightJustify(cfg.format(cfg.fixed(std::log(1ll << i))));
    }
    file << "};\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     Variables\n";
    file << "    // ----============================================================----\n\n";

    file << "    // Real and imaginary part of L_k.\n"
         << "    " << cfg.dataType() << " rLk_r, rLk_i;\n";
    file << "    // Real and imaginary part of E_k.\n"
         << "    " << cfg.dataType() << " rEk_r, rEk_i;\n";
    file << "\n";

    file << "    // Fractions of L_k. Needed only for reduction.\n"
         << "    " << cfg.dataType() << " rLk_r_div_2, rLk_r_div_4, rLk_i_div_2;\n\n";

    file << "    bit [4:0] r5Rotated;\n"
         << "    bit [" << reductionShiftCounter.numBits - 1 << ":0] " << reductionShiftCounter.name << ";\n"
         << "    bit rReductionScalingByMultiplication;\n\n";

    file << "    wire `AssignDelay wL_r_underflow = rLk_r["
         << cfg.HB - 1 << ":" << -cfg.FB << "] < CONST_UNDERFLOW_THRESHOLD;\n";
    file << "    wire `AssignDelay wL_r_overflow  = ";
    if (overflowThreshold > 0)
    {
        file << "rLk_r[" << cfg.HB - 1 << ":" << -cfg.FB << "] > CONST_OVERFLOW_THRESHOLD;\n\n";
    }
    else
    {
        file << "1'b0;\n\n";
    }

    file << "    // Signals if real(L) is below the convergence bound. Note that real(L)> = 0.\n"
         << "    wire `AssignDelay wL_r_too_low  = rLk_r[" << cfg.HB - 1 << ":" << -cfg.FB << "] < CONST_0_64;\n"
         << "    // Signals if real(L) is above the convergence bound. Note that real(L)> = 0.\n"
         << "    wire `AssignDelay wL_r_too_high = rLk_r[" << cfg.HB - 1 << ":" << -cfg.FB << "] > CONST_1_40;\n\n";

    file << "    wire `AssignDelay wE_i_too_low  = rEk_i <= -CONST_PI;\n"
         << "    wire `AssignDelay wE_i_too_high = rEk_i > CONST_PI;\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     State machine\n";
    file << "    // ----============================================================----\n\n";

    file << "    typedef enum\n"
         << "    {\n"
         << "        State_WaitForInput,\n"

         << "        State_Reduce_RotateToFirstQuadrant,\n"
         << "        State_Reduce_RotateBelowDiagonal,\n"
         << "        State_Reduce_ComputeFractions,\n"
         << "        State_Reduce_RotateToConvergenceAngle,\n"
         << "        State_Reduce_ScaleMagnitude,\n"

         << "        State_Compute_InitStage1,\n"
         << "        State_Compute_InitStage2,\n"
         << "        State_Compute_Stage1,\n"
         << "        State_Compute_Stage2,\n"
         << "        State_Compute_Stage3,\n"
         << "        State_Compute_Stage4,\n"

         << "        State_Expand_FirstQuadrantRotation,\n"
         << "        State_Expand_RotationByMinusPiDiv4_1,\n"
         << "        State_Expand_RotationByMinusPiDiv4_2,\n"
         << "        State_Expand_RotationByAtanMinus1Div2,\n"
         << "        State_Expand_RescaleAndNormalizeAngle,\n"

         << "        State_SetOutOfRange,\n"
         << "        State_Output\n"
         << "    } State_t;\n\n";
    file << "    State_t state, nextState;\n\n";

    file << "    bit [" << iterationCounter.numBits - 1 << ":0] " << iterationCounter.name << ";\n";
    file << "    bit [" << lookupOffset.numBits - 1 << ":0] " << lookupOffset.name << ";\n\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "        begin\n"
         << "            state <= `AssignDelay State_WaitForInput;\n"
         << "        end\n"
         << "        else\n"
         << "        begin\n"
         << "            state <= `AssignDelay nextState;\n"
         << "        end\n"
         << "    end\n\n";

    file << "    // Next state logic.\n";
    file << "    always_comb\n"
         << "    begin\n"
         << "        unique case (state)\n"
         << "            State_WaitForInput:\n"
         << "                nextState <= inputValid ? State_Reduce_RotateToFirstQuadrant : State_WaitForInput;\n\n"

         << "            State_Reduce_RotateToFirstQuadrant:\n"
         << "                nextState <= State_Reduce_RotateBelowDiagonal;\n"
         << "            State_Reduce_RotateBelowDiagonal:\n"
         << "                nextState <= State_Reduce_ComputeFractions;\n"
         << "            State_Reduce_ComputeFractions:\n"
         << "                nextState <= State_Reduce_RotateToConvergenceAngle;\n"
         << "            State_Reduce_RotateToConvergenceAngle:\n"
         << "                nextState <= State_Reduce_ScaleMagnitude;\n"
         << "            State_Reduce_ScaleMagnitude:\n"
         << "            begin\n"
         << "                if (wL_r_underflow || wL_r_overflow)\n"
         << "                    nextState <= State_SetOutOfRange;\n"
         << "                else if (wL_r_too_low || wL_r_too_high)\n"
         << "                    nextState <= State_Reduce_ScaleMagnitude;\n"
         << "                else\n"
         << "                    nextState <= State_Compute_InitStage1;\n"
         << "            end\n\n"

         << "            State_Compute_InitStage1:\n"
         << "                nextState <= State_Compute_InitStage2;\n"
         << "            State_Compute_InitStage2:\n"
         << "                nextState <= State_Compute_Stage1;\n"
         << "            State_Compute_Stage1:\n"
         << "                nextState <= State_Compute_Stage2;\n"
         << "            State_Compute_Stage2:\n"
         << "                nextState <= State_Compute_Stage3;\n"
         << "            State_Compute_Stage3:\n"
         << "                nextState <= State_Compute_Stage4;\n"
         << "            State_Compute_Stage4:\n"
         << "                nextState <= " << iterationCounter.name << " == " << iterationCounter.numBits << "'d" << numIter << " ? State_Expand_FirstQuadrantRotation : State_Compute_Stage1;\n\n"

         << "            State_Expand_FirstQuadrantRotation:\n"
         << "                nextState <= State_Expand_RotationByMinusPiDiv4_1;\n"
         << "            State_Expand_RotationByMinusPiDiv4_1:\n"
         << "                nextState <= State_Expand_RotationByMinusPiDiv4_2;\n"
         << "            State_Expand_RotationByMinusPiDiv4_2:\n"
         << "                nextState <= State_Expand_RotationByAtanMinus1Div2;\n"
         << "            State_Expand_RotationByAtanMinus1Div2:\n"
         << "                nextState <= State_Expand_RescaleAndNormalizeAngle;\n"
         << "            State_Expand_RescaleAndNormalizeAngle:\n"
         << "                nextState <= (wE_i_too_low || wE_i_too_high) ? State_Expand_RescaleAndNormalizeAngle : State_Output;\n\n"

         << "            State_SetOutOfRange:\n"
         << "                nextState <= State_Output;\n\n"

         << "            State_Output:\n"
         << "                nextState <= State_WaitForInput;\n\n"

         << "            default:\n"
         << "                nextState <= State_WaitForInput;\n"
         << "        endcase\n"
         << "    end\n\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "        begin\n"
         << "            " << iterationCounter.name << " <= `AssignDelay '0;\n"
         << "            " << lookupOffset.name << " <= `AssignDelay '0;\n"
         << "        end\n"
         << "        else if (state == State_Compute_InitStage1)\n"
         << "        begin\n"
         << "            " << iterationCounter.name << " <= `AssignDelay '0;\n"
         << "            " << lookupOffset.name << " <= `AssignDelay '0;\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage3)\n"
         << "        begin\n"
         << "            " << iterationCounter.name << " <= `AssignDelay " << iterationCounter.name << " + 1'd1;\n"
         << "            " << lookupOffset.name << " <= `AssignDelay " << lookupOffset.name << " + 4'd8;\n"
         << "        end\n"
         << "    end\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     d_k\n";
    file << "    // ----============================================================----\n\n";

    file << "    typedef enum bit [2:0]\n";
    file << "    {\n";
    file << "        Dk_N = 3'b001,\n";
    file << "        Dk_0 = 3'b010,\n";
    file << "        Dk_P = 3'b100\n";
    file << "    } DkValue_t;\n\n";

    file << "    DkValue_t rD_r, rD_i;\n";

    const int tSize = cfg.IB + 4;

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "        begin\n"
         << "            rD_r <= `AssignDelay Dk_0;\n"
         << "            rD_i <= `AssignDelay Dk_0;\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage1)\n"
         << "        begin\n"
         << "            if (" << iterationCounter.name << " == '0)\n"
         << "            begin\n"
         << "                // Choose d_1 = 0.\n"
         << "                rD_r <= `AssignDelay Dk_0;\n"
         << "                rD_i <= `AssignDelay Dk_0;\n"
         << "            end\n"
         << "            else\n"
         << "            begin\n"
         << "                // Determine the real part of d_k from floor(L_k^r * 16).\n"
         << "                if (signed'(rLk_r[" << cfg.HB << ":-4]) <= -" << tSize << "'sd8)\n"
         << "                    rD_r <= `AssignDelay Dk_P;\n"
         << "                else if (signed'(rLk_r[" << cfg.HB << ":-4]) >= " << tSize << "'sd8)\n"
         << "                    rD_r <= `AssignDelay Dk_N;\n"
         << "                else\n"
         << "                    rD_r <= `AssignDelay Dk_0;\n"

         << "                // Determine the imaginary part of d_k from floor(L_k^i * 16).\n"
         << "                if (signed'(rLk_i[" << cfg.HB << ":-4]) <= -" << tSize << "'sd8)\n"
         << "                    rD_i <= `AssignDelay Dk_P;\n"
         << "                else if (signed'(rLk_i[" << cfg.HB << ":-4]) >= " << tSize << "'sd8)\n"
         << "                    rD_i <= `AssignDelay Dk_N;\n"
         << "                else\n"
         << "                    rD_i <= `AssignDelay Dk_0;\n"
         << "            end\n"
         << "        end\n"
         << "    end\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     Lookup\n";
    file << "    // ----============================================================----\n\n";

    file << "    " << cfg.dataType() << " rLnValue, rArctanValue;\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "        begin\n"
         << "            rLnValue     <= `AssignDelay '0;\n"
         << "            rArctanValue <= `AssignDelay '0;\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage2)\n"
         << "        begin\n"
         << "            // Note: We should fetch rLnValue in this stage.\n"
         << "            // We fetch into rArctanValue and copy to rLnValue later on in\n"
         << "            // the hope that during synthesis a block ram will be inferred.\n"
         << "            if (rD_r & Dk_N)\n"
         << "            begin\n"
         << "                if (rD_i & Dk_0) rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd0];\n"
         << "                else             rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd1];\n"
         << "            end\n"
         << "            else if (rD_r & Dk_0)\n"
         << "            begin\n"
         << "                if (rD_i & Dk_0) rArctanValue <= `AssignDelay '0;\n"
         << "                else             rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd2];\n"
         << "            end\n"
         << "            else\n"
         << "            begin\n"
         << "                if (rD_i & Dk_0) rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd3];\n"
         << "                else             rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd4];\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage3)\n"
         << "        begin\n"
         << "            rLnValue <= `AssignDelay rArctanValue;\n"
         << "            if      (rD_r & Dk_N) rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd5];\n"
         << "            else if (rD_r & Dk_0) rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd6];\n"
         << "            else                  rArctanValue <= `AssignDelay LUT_BKM[" << lookupOffset.name << " + 4'd7];\n"
         << "        end\n"
         << "    end\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     H_k\n";
    file << "    // ----============================================================----\n\n";

    file << "    " << cfg.dataType() << " rH_r, rH_i;\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "        begin\n"
         << "            rH_r <= `AssignDelay '0;\n"
         << "            rH_i <= `AssignDelay '0;\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage2)\n"
         << "        begin\n"
         << "            // H_k^r = (L_k^r * d^r - L_k^i * d^i)\n"
         << "            // H_k^i = (L_k^r * d^i + L_k^i * d^r)\n"
         << "            if (rD_r & Dk_N)\n"
         << "            begin\n"
         << "                if      (rD_i & Dk_N) rH_r <= `AssignDelay  rLk_i - rLk_r;\n"
         << "                else if (rD_i & Dk_0) rH_r <= `AssignDelay -rLk_r;\n"
         << "                else                  rH_r <= `AssignDelay -rLk_r - rLk_i;\n"

         << "                if      (rD_i & Dk_N) rH_i <= `AssignDelay -rLk_r - rLk_i;\n"
         << "                else if (rD_i & Dk_0) rH_i <= `AssignDelay -rLk_i;\n"
         << "                else                  rH_i <= `AssignDelay  rLk_r - rLk_i;\n"
         << "            end\n"
         << "            else if (rD_r & Dk_0)\n"
         << "            begin\n"
         << "                if      (rD_i & Dk_N) rH_r <= `AssignDelay  rLk_i;\n"
         << "                else if (rD_i & Dk_0) rH_r <= `AssignDelay  '0;\n"
         << "                else                  rH_r <= `AssignDelay -rLk_i;\n"

         << "                if      (rD_i & Dk_N) rH_i <= `AssignDelay -rLk_r;\n"
         << "                else if (rD_i & Dk_0) rH_i <= `AssignDelay  '0;\n"
         << "                else                  rH_i <= `AssignDelay  rLk_r;\n"
         << "            end\n"
         << "            else\n"
         << "            begin\n"
         << "                if      (rD_i & Dk_N) rH_r <= `AssignDelay  rLk_r + rLk_i;\n"
         << "                else if (rD_i & Dk_0) rH_r <= `AssignDelay  rLk_r;\n"
         << "                else                  rH_r <= `AssignDelay  rLk_r - rLk_i;\n"

         << "                if      (rD_i & Dk_N) rH_i <= `AssignDelay  rLk_i - rLk_r;\n"
         << "                else if (rD_i & Dk_0) rH_i <= `AssignDelay  rLk_i;\n"
         << "                else                  rH_i <= `AssignDelay  rLk_r + rLk_i;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage3)\n"
         << "        begin\n"
         << "            // H_k = H_k * 2^(-k+1)\n"
         << "            rH_r <= `AssignDelay rH_r >>> " << iterationCounter.name << ";\n"
         << "            rH_i <= `AssignDelay rH_i >>> " << iterationCounter.name << ";\n"
         << "        end\n"
         << "        else if (state == State_Expand_RotationByAtanMinus1Div2)\n"
         << "        begin\n"
         << "            rH_r <= `AssignDelay LUT_LOG_POW2[" << reductionShiftCounter.name << "];\n"
         << "        end\n"
         << "        else if (state == State_Expand_RescaleAndNormalizeAngle)\n"
         << "        begin\n"
         << "            rH_r <= `AssignDelay '0;\n"
         << "        end\n"
         << "    end\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     E_k\n";
    file << "    // ----============================================================----\n\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "        begin\n"
         << "            rEk_r <= `AssignDelay '0;\n"
         << "            rEk_i <= `AssignDelay '0;\n"
         << "        end\n"
         << "        else if (state == State_WaitForInput)\n"
         << "        begin\n"
         << "            if (inputValid)\n"
         << "            begin\n"
         << "                rEk_r <= `AssignDelay E_in_r;\n"
         << "                rEk_i <= `AssignDelay E_in_i;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage4)\n"
         << "        begin\n"
         << "            // E_k = E_k - Delta E_k\n"
         << "            rEk_r <= `AssignDelay rEk_r - rLnValue;\n"
         << "            if      (rD_i & Dk_N) rEk_i <= `AssignDelay rEk_i + rArctanValue;\n"
         << "            else if (rD_i & Dk_P) rEk_i <= `AssignDelay rEk_i - rArctanValue;\n"
         << "        end\n"

         << "        else if (state == State_Expand_FirstQuadrantRotation)\n"
         << "        begin\n"
         << "            if (r5Rotated[1:0] == 2'b11)\n"
         << "                rEk_i <= `AssignDelay rEk_i - CONST_PI;\n"
         << "            else if (r5Rotated[1:0] == 2'b01)\n"
         << "                rEk_i <= `AssignDelay rEk_i + CONST_PI_2;\n"
         << "            else if (r5Rotated[1:0] == 2'b10)\n"
         << "                rEk_i <= `AssignDelay rEk_i - CONST_PI_2;\n"
         << "        end\n"
         << "        else if (state == State_Expand_RotationByMinusPiDiv4_1)\n"
         << "        begin\n"
         << "            if (r5Rotated[2])\n"
         << "            begin\n"
         << "                rEk_r <= `AssignDelay rEk_r + CONST_LOG_1_SQRT2;\n"
         << "                rEk_i <= `AssignDelay rEk_i + CONST_PI_4;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Expand_RotationByMinusPiDiv4_2)\n"
         << "        begin\n"
         << "            if (r5Rotated[3])\n"
         << "            begin\n"
         << "                rEk_r <= `AssignDelay rEk_r + CONST_LOG_1_SQRT2;\n"
         << "                rEk_i <= `AssignDelay rEk_i + CONST_PI_4;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Expand_RotationByAtanMinus1Div2)\n"
         << "        begin\n"
         << "            if (r5Rotated[4])\n"
         << "            begin\n"
         << "                rEk_r <= `AssignDelay rEk_r + CONST_LOG_2_SQRT5;\n"
         << "                rEk_i <= `AssignDelay rEk_i + CONST_ATAN_1_2;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Expand_RescaleAndNormalizeAngle)\n"
         << "        begin\n"
         << "            if (rReductionScalingByMultiplication)\n"
         << "                rEk_r <= `AssignDelay rEk_r - rH_r;\n"
         << "            else\n"
         << "                rEk_r <= `AssignDelay rEk_r + rH_r;\n"

         << "            if (wE_i_too_low)\n"
         << "                rEk_i <= `AssignDelay rEk_i + CONST_PI;\n"
         << "            else if (wE_i_too_high)\n"
         << "                rEk_i <= `AssignDelay rEk_i - CONST_PI;\n"
         << "        end\n"

         << "        else if (state == State_SetOutOfRange)\n"
         << "        begin\n"
         << "            if (wL_r_underflow)\n"
         << "                rEk_r <= `AssignDelay {1'b1, '0};\n"
         << "            else\n"
         << "                rEk_r <= `AssignDelay {1'b0, '1};\n"
         << "            rEk_i <= `AssignDelay '0;\n"
         << "        end\n"

         << "    end\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     L_k\n";
    file << "    // ----============================================================----\n\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "        begin\n"
         << "            rLk_r <= `AssignDelay '0;\n"
         << "            rLk_i <= `AssignDelay '0;\n"
         << "            r5Rotated <= `AssignDelay '0;\n"
         << "            " << reductionShiftCounter.name << " <= `AssignDelay '0;\n"
         << "        end\n"
         << "        else if (state == State_WaitForInput)\n"
         << "        begin\n"
         << "            if (inputValid)\n"
         << "            begin\n"
         << "                rLk_r <= `AssignDelay L_in_r;\n"
         << "                rLk_i <= `AssignDelay L_in_i;\n"
         << "                r5Rotated <= `AssignDelay '0;\n"
         << "                " << reductionShiftCounter.name << " <= `AssignDelay '0;\n"
         << "            end\n"
         << "        end\n"

         << "        else if (state == State_Reduce_RotateToFirstQuadrant)\n"
         << "        begin\n"
         << "            if (rLk_r < 0 && rLk_i < 0)\n"
         << "            begin\n"
         << "                // Rotate by pi.\n"
         << "                rLk_r <= `AssignDelay -rLk_r;\n"
         << "                rLk_i <= `AssignDelay -rLk_i;\n"
         << "                r5Rotated[1:0] <= `AssignDelay 2'b11;\n"
         << "            end\n"
         << "            else if (rLk_r < 0)\n"
         << "            begin\n"
         << "                // Rotate by -pi/2.\n"
         << "                rLk_r <= `AssignDelay  rLk_i;\n"
         << "                rLk_i <= `AssignDelay -rLk_r;\n"
         << "                r5Rotated[1:0] <= `AssignDelay 2'b01;\n"
         << "            end\n"
         << "            else if (rLk_i < 0)\n"
         << "            begin\n"
         << "                // Rotate by pi/2.\n"
         << "                rLk_r <= `AssignDelay -rLk_i;\n"
         << "                rLk_i <= `AssignDelay  rLk_r;\n"
         << "                r5Rotated[1:0] <= `AssignDelay 2'b10;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Reduce_RotateBelowDiagonal)\n"
         << "        begin\n"
         << "            if (rLk_i > rLk_r)\n"
         << "            begin\n"
         << "                // Rotate by -pi/4. The scaling by 1/sqrt(2) is omitted.\n"
         << "                rLk_r <= `AssignDelay rLk_r + rLk_i;\n"
         << "                rLk_i <= `AssignDelay rLk_i - rLk_r;\n"
         << "                r5Rotated[2] <= `AssignDelay 1'b1;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Reduce_ComputeFractions)\n"
         << "        begin\n"
         << "            rLk_r_div_2 <= `AssignDelay rLk_r >>> 1;\n"
         << "            rLk_r_div_4 <= `AssignDelay rLk_r >>> 2;\n"
         << "            rLk_i_div_2 <= `AssignDelay rLk_i >>> 1;\n"
         << "        end\n"
         << "        else if (state == State_Reduce_RotateToConvergenceAngle)\n"
         << "        begin\n"
         << "            if (rLk_i > rLk_r_div_2)\n"
         << "            begin\n"
         << "                // Rotate by -pi/4. The scaling by 1/sqrt(2) is omitted.\n"
         << "                rLk_r <= `AssignDelay rLk_r + rLk_i;\n"
         << "                rLk_i <= `AssignDelay rLk_i - rLk_r;\n"
         << "                r5Rotated[3] <= `AssignDelay 1'b1;\n"
         << "            end\n"
         << "            else if (rLk_i > rLk_r_div_4)\n"
         << "            begin\n"
         << "                // Rotate by atan(-1/2). The scaling by 2/sqrt(5) is omitted.\n"
         << "                rLk_r <= `AssignDelay rLk_r + rLk_i_div_2;\n"
         << "                rLk_i <= `AssignDelay rLk_i - rLk_r_div_2;\n"
         << "                r5Rotated[4] <= `AssignDelay 1'b1;\n"
         << "            end\n"
         << "        end\n"
         << "        else if (state == State_Reduce_ScaleMagnitude)\n"
         << "        begin\n"
         << "            if (wL_r_too_low)\n"
         << "            begin\n"
         << "                rReductionScalingByMultiplication <= `AssignDelay 1'b1;\n"
         << "                rLk_r <= `AssignDelay rLk_r <<< 1'd1;\n"
         << "                rLk_i <= `AssignDelay rLk_i <<< 1'd1;\n"
         << "                " << reductionShiftCounter.name << " <= `AssignDelay " << reductionShiftCounter.name << " + 1'd1;\n"
         << "            end\n"
         << "            else if (wL_r_too_high)\n"
         << "            begin\n"
         << "                rReductionScalingByMultiplication <= `AssignDelay 1'b0;\n"
         << "                rLk_r <= `AssignDelay rLk_r >>> 1'd1;\n"
         << "                rLk_i <= `AssignDelay rLk_i >>> 1'd1;\n"
         << "                " << reductionShiftCounter.name << " <= `AssignDelay " << reductionShiftCounter.name << " + 1'd1;\n"
         << "            end\n"
         << "        end\n"

         << "        else if (state == State_Compute_InitStage1)\n"
         << "        begin\n"
         << "            rLk_r <= `AssignDelay rLk_r - " << cfg.format(cfg.fixed(1)) << ";\n"
         << "        end\n"
         << "        else if (state == State_Compute_InitStage2)\n"
         << "        begin\n"
         << "            rLk_r <= `AssignDelay rLk_r <<< 1'd1;\n"
         << "            rLk_i <= `AssignDelay rLk_i <<< 1'd1;\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage2)\n"
         << "        begin\n"
         << "            // L_k = L_k + d_k\n"
         << "            if      (rD_r & Dk_N) rLk_r <= `AssignDelay rLk_r - " << cfg.format(cfg.fixed(1)) << ";\n"
         << "            else if (rD_r & Dk_P) rLk_r <= `AssignDelay rLk_r + " << cfg.format(cfg.fixed(1)) << ";\n\n"
         << "            if      (rD_i & Dk_N) rLk_i <= `AssignDelay rLk_i - " << cfg.format(cfg.fixed(1)) << ";\n"
         << "            else if (rD_i & Dk_P) rLk_i <= `AssignDelay rLk_i + " << cfg.format(cfg.fixed(1)) << ";\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage3)\n"
         << "        begin\n"
         << "            // L_k = 2 * L_k\n"
         << "            rLk_r <= `AssignDelay rLk_r <<< 1'd1;\n"
         << "            rLk_i <= `AssignDelay rLk_i <<< 1'd1;\n"
         << "        end\n"
         << "        else if (state == State_Compute_Stage4)\n"
         << "        begin\n"
         << "            // L_k = L_k + H_k\n"
         << "            rLk_r <= `AssignDelay rLk_r + rH_r;\n"
         << "            rLk_i <= `AssignDelay rLk_i + rH_i;\n"
         << "        end\n"
         << "    end\n\n";

    file << "    // ----============================================================----\n";
    file << "    //     Output\n";
    file << "    // ----============================================================----\n\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "            outputValid <= `AssignDelay 1'b0;\n"
         << "        else if (state == State_Output)\n"
         << "            outputValid <= `AssignDelay 1'b1;\n"
         << "        else\n"
         << "            outputValid <= `AssignDelay 1'b0;\n"
         << "    end\n\n";

    file << "    always_ff @(posedge clock or negedge reset_n)\n"
         << "    begin\n"
         << "        if (!reset_n)\n"
         << "            outOfRange <= `AssignDelay 1'b0;\n"
         << "        else if (state == State_WaitForInput)\n"
         << "            outputValid <= `AssignDelay 1'b0;\n"
         << "        else if (state == State_SetOutOfRange)\n"
         << "            outOfRange <= `AssignDelay 1'b1;\n"
         << "    end\n\n";

    file << "    assign E_out_r = rEk_r;\n"
         << "    assign E_out_i = rEk_i;\n"
         << "    assign L_out_r = rLk_r;\n"
         << "    assign L_out_i = rLk_i;\n\n";

    file << "endmodule\n\n";
}

void testbench(const Config& cfg)
{
    std::ofstream file("testbench.sv");
    file.precision(numeric_limits<double>::digits10 + 2);

    string pow;
    {
        stringstream ss;
        ss << "2**";
        ss << cfg.FB;
        pow = ss.str();
    }

    file << "`timescale 1ns/1ns\n\n";

    file << "module top;\n\n";

    file << "bit clock;\n"
         << "bit reset_n;\n"
         << "bit signed [3:-12] L_in_r;\n"
         << "bit signed [3:-12] L_in_i;\n"
         << "bit signed [3:-12] E_in_r;\n"
         << "bit signed [3:-12] E_in_i;\n"
         << "bit inputValid;\n"
         << "bit signed [3:-12] L_out_r;\n"
         << "bit signed [3:-12] L_out_i;\n"
         << "bit signed [3:-12] E_out_r;\n"
         << "bit signed [3:-12] E_out_i;\n"
         << "bit outputValid;\n"
         << "bit outOfRange;\n\n";

    file << "BkmLComputer DUT(.*);\n\n";

    file << "// System clock\n"
         << "initial\n"
         << "begin\n"
         << "    clock = 0;\n"
         << "    forever\n"
         << "        #500ns clock = ~clock;\n"
         << "end\n\n";

    file << "event dutFinished;\n\n";

    file << "real expectedMagnitude;\n"
         << "real expectedPhase;\n\n";

    file << "task runOnce(real magnitude, real phase);\n"
         << "begin\n"
         << "    while (phase <= " << -pi << ")\n"
         << "    begin\n"
         << "        phase += " << 2 * pi << ";\n"
         << "    end\n"
         << "    while (phase > " << pi << ")\n"
         << "    begin\n"
         << "        phase -= " << 2 * pi << ";\n"
         << "    end\n\n"

         << "    @(negedge clock);\n"
         << "    expectedMagnitude = magnitude;\n"
         << "    expectedPhase = phase;\n"
         << "    L_in_r <= $rtoi(expectedMagnitude * $cos(expectedPhase) * " << pow << ");\n"
         << "    L_in_i <= $rtoi(expectedMagnitude * $sin(expectedPhase) * " << pow << ");\n"
         << "    E_in_r <= 0;\n"
         << "    E_in_i <= 0;\n"
         << "    inputValid <= 1;\n"
         << "    @(negedge clock);\n"
         << "    inputValid <= 0;\n"
         << "    @dutFinished;\n"
         << "    @(negedge clock);\n"
         << "end\n"
         << "endtask\n\n";

    file << "task runCircle(real magnitude);\n"
         << "begin\n"
         << "    automatic real pi = " << pi << ";\n"
         << "    for (int idx = -32 + 1; idx <= 32; ++idx)\n"
         << "        runOnce(magnitude, pi / 32 * idx);\n"
         << "end\n"
         << "endtask\n\n";

    file << "initial\n"
         << "begin\n"
         << "    reset_n = 0;\n"
         << "    #1.5us;\n"
         << "    reset_n = 1;\n"
         << "    #2us;\n"
         << "    @(negedge clock);\n\n";

    file << "    // Check unit circle\n"
         << "    runCircle(1);\n\n";

    file << "    // Check around convergence bounds\n"
         << "    for (int diff = -10; diff <= 10; ++diff)\n"
         << "        runCircle(0.64 + $itor(diff) / " << pow << ");\n"
         << "    for (int diff = -10; diff <= 10; ++diff)\n"
         << "        runCircle(1.4 + $itor(diff) / " << pow << ");\n\n";

    file << "    $display(\"==== Done ====\");\n";

    file << "end\n\n";

    file << "initial\n"
         << "begin\n"
         << "    automatic real log;\n"
         << "    automatic real logDiff;\n"
         << "    automatic real phi;\n"
         << "    automatic real phiDiff;\n\n";

    file << "    $display(\"|--LOGARITHM-----------------------------------------------|          "
                           "|--PHASE---------------------------------------------------|\");\n";
    file << "    $display(\"%10s  %10s  %10s  %10s  %10s  %20s  %10s  %10s  %10s  %10s\",\n"
            "             \"fixed\", \"actual\", \"expected\", \"delta\", \"fix-delta\",\n"
            "             \"fixed\", \"actual\", \"expected\", \"delta\", \"fix-delta\");\n";

    file << "    forever\n"
         << "    begin\n"
         << "        @(posedge outputValid);\n"
         << "        log = $itor(E_out_r) / " << pow << ";\n"
         << "        logDiff = log - $ln(expectedMagnitude);\n\n"
         << "        phi = $itor(E_out_i) / " << pow << ";\n"
         << "        phiDiff = phi - expectedPhase;\n\n"

         //<< "$display("Finished at %t", $time);
         << "        $display(\"%10d, %10f, %10f, %10f, %10f, %20d, %10f, %10f, %10f, %10f\",\n"
            "                 E_out_r, log, $ln(expectedMagnitude), logDiff, logDiff * " << pow << ",\n"
            "                 E_out_i, phi, expectedPhase, phiDiff, phiDiff * " << pow << ");\n\n"

         << "        ->dutFinished;\n"

         << "    end\n"
         << "end\n\n";

    file << "endmodule\n\n";
}

int main()
{
    std::ofstream file("bkmcomputer.sv");

    ConfigBuilder builder;
    builder.setNumBits(16);
    builder.setNumFractionBits(12);
    builder.setUnderflowThreshold(2);
    Config cfg = builder.create();
    lModeComputer(file, cfg);
    testbench(cfg);
}
