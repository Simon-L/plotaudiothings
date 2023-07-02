#include <iostream>
#include <string>
#include <cstdio>

#include "sst/basic-blocks/modulators/ADSREnvelope.h"
#include "digital.hpp"
#include "CParamSmooth.hpp"
#include <sciplot/sciplot.hpp>
using namespace sciplot;

std::vector<float> out;
std::vector<float> out2;
std::vector<float> out3;
std::vector<float> time_x;


static constexpr int tbs{16};
struct SampleSRProvider
{
    double samplerate{44100}, sampleRateInv{1.f / samplerate};
    float envelope_rate_linear_nowrap(float f) const { return tbs * sampleRateInv * pow(2.f, -f); }
} srp;


rack::dsp::PulseGenerator trigger;

auto main() -> int
{
    auto adsr = sst::basic_blocks::modulators::ADSREnvelope<SampleSRProvider, tbs>(&srp);
    
    float _A = 0.001;
    float _D = 0.75;
    float _S = 0.0f;
    float _R = 4 * 0.10;
    
    float deltaTime = 1./44100.;
    float length = 0.3;
    
    trigger.trigger(0.650);
    adsr.attackFrom(0.0, _A, 1, false); // initial, attacktime, ashp, digital?
    
    CParamSmooth pasm(20, 44100);
    
    for (size_t s = 0; s < 44100; s++) {
        auto b = trigger.process(deltaTime);
            
        adsr.process(_A, _D, _S, _R, 1, 1, 1, b); // a, d, s, r, ashp, dshp, rshp, gateActive
        out.push_back(adsr.output);
        out3.push_back(pasm.process(adsr.output));
        out2.push_back(s < (0.16 * 44100) ? 1.0 : 0.0);
        time_x.push_back(s * deltaTime);
        // if (s < 100) printf("%f\n", adsr.output);
        // out[s] = b ? 1.0 : -1.0;
    }
    
    Plot2D plot;
    plot.xlabel("s");
    plot.ylabel("env");
    plot.xrange(0.0, 44100 * deltaTime);
    plot.yrange(0.0, 1.0);
    plot.drawCurve(time_x, out);
    // plot.drawCurve(time_x, out2);
    plot.drawCurve(time_x, out3);
    Figure fig = {{plot}};
    Canvas canvas = {{fig}};
    canvas.size(1000,1000);
    canvas.show();
    
    return 0;
}
