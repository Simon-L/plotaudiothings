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
#include "values_reso.hpp"
#include "Resonator.hpp"

auto main(int argc, const char *argv[]) -> int
{   
    
    const juce::File optionsFileJ = juce::File(argv[1]);
    chowdsp::GenericTweaksFile<false> options;
    options.initialise (optionsFileJ, 1);
    
    auto freq = options.getProperty<float> ("freq");
    auto q = options.getProperty<float> ("q");
    auto gain = options.getProperty<float> ("gain");
    auto xmin = options.getProperty<float> ("xmin");
    auto xmax = options.getProperty<float> ("xmax");
    auto ymin = options.getProperty<float> ("ymin");
    auto ymax = options.getProperty<float> ("ymax");
    auto Rfb = options.getProperty<float> ("Rfb");
    auto R_g = options.getProperty<float> ("R_g");
    auto C = options.getProperty<float> ("C");
    
    printf("Options: %f %f %f\n", freq, q, gain);
    
    chowdsp::ButterworthFilter< 3, chowdsp::ButterworthFilterType::Highpass, float> fi;
    HatResonatorWDF reso;
    reso.prepare(48000);
    // reso.setParameters(Rfb, R_g, C);
    chowdsp::Gain<float> preGain;
    preGain.setGainLinear(gain);
    
    
    chowdsp::SpectrumPlotBase base {
            chowdsp::SpectrumPlotParams {
                500.0f,
                20000.0f,
                -30.0f,
                30.0f }
    };
    chowdsp::GenericFilterPlotter plotter { base, {} };
    plotter.runFilterCallback = [&fi, &preGain, gain, freq, q, &reso] (const float* in, float* out, int N)
    {   
        std::copy (in, in + N, out);
        preGain.prepare({48000.0, N, 1});
        preGain.process (chowdsp::BufferView { out, N });
        // fi.reset();
        reso.reset();
        reso.prepare(48000);
        // reso.setParameters(8000.0, 0.06);
        // fi.prepare(1);
        // fi.calcCoefs(freq, q, 48000.0);
        // fi.processBlock (chowdsp::BufferView { out, N });
        for (size_t i = 0; i < N; i++) {
            out[i] = reso.processSample(out[i]);
        }
    };

    const auto [freqAxis, magAxis] = plotter.plotFilterMagnitudeResponse();    
    
    Plot2D plot;
    plot.xlabel("Frequency (Hz)");
    plot.ylabel("Magnitude (dB)");
    plot.drawCurve(freqAxis, magAxis).label("chowdsp::ButterworthFilter");
    plot.drawCurve(freqValues, dbValues).label("LTspice sim");
    plot.xtics().logscale();
    plot.xrange(xmin, xmax);
    plot.yrange(ymin, ymax);
    Figure fig = {{plot}};
    Canvas canvas = {{fig}};
    canvas.size(1300,1300);
    canvas.show();
    
    return 0;
}
