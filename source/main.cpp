#include <iostream>
#include <string>
#include <cstdio>

#include "sst/basic-blocks/modulators/ADSREnvelope.h"
#include "digital.hpp"
#include <sciplot/sciplot.hpp>
using namespace sciplot;

#include "chowdsp_dsp_utils/chowdsp_dsp_utils.h"
#include "chowdsp_sources/chowdsp_sources.h"
#include "chowdsp_visualizers/chowdsp_visualizers.h"
#include <chowdsp_plugin_utils/chowdsp_plugin_utils.h>

// #include "values_ch.hpp"
#include "values_oh.hpp"

auto main(int argc, const char *argv[]) -> int
{   
    
    const juce::File optionsFileJ = juce::File(argv[1]);
    chowdsp::GenericTweaksFile<false> options;
    options.initialise (optionsFileJ, 1);
    
    auto freq = options.getProperty<float> ("freq");
    auto q = options.getProperty<float> ("q");
    auto gain = options.getProperty<float> ("gain");
    
    printf("Options: %f %f %f\n", freq, q, gain);
    
    chowdsp::ButterworthFilter< 2, chowdsp::ButterworthFilterType::Highpass, float> fi;
    chowdsp::Gain<float> preGain;
    preGain.setGainLinear(gain);
    
    fi.prepare(1);
    fi.calcCoefs(freq, q, 48000);
    
    chowdsp::SpectrumPlotBase base {
            chowdsp::SpectrumPlotParams {
                500.0f,
                20000.0f,
                -30.0f,
                30.0f }
    };
    chowdsp::GenericFilterPlotter plotter { base, {} };
    plotter.runFilterCallback = [&fi, &preGain, gain] (const float* in, float* out, int N)
    {   
        std::copy (in, in + N, out);
        preGain.prepare({48000, N, 1});
        preGain.process (chowdsp::BufferView { out, N });
        fi.processBlock (chowdsp::BufferView { out, N });
    };

    const auto [freqAxis, magAxis] = plotter.plotFilterMagnitudeResponse();    
    
    Plot2D plot;
    plot.xlabel("Frequency (Hz)");
    plot.ylabel("Magnitude (dB)");
    plot.xrange(10000.0f, 15000.0f);
    plot.yrange(-30.0f, 8.0f);
    plot.drawCurve(freqAxis, magAxis).label("chowdsp::ButterworthFilter");
    plot.drawCurve(freqValues, dbValues).label("LTspice sim");
    plot.xtics().logscale();
    Figure fig = {{plot}};
    Canvas canvas = {{fig}};
    canvas.size(1000,1000);
    canvas.show();
    
    return 0;
}
