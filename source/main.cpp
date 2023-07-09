#include <iostream>
#include <string>
#include <cstdio>
#include <regex>

#include "sst/basic-blocks/modulators/ADSREnvelope.h"
#include "digital.hpp"
#include <sciplot/sciplot.hpp>
using namespace sciplot;

#include "chowdsp_dsp_utils/chowdsp_dsp_utils.h"
#include "chowdsp_sources/chowdsp_sources.h"
#include "chowdsp_visualizers/chowdsp_visualizers.h"
#include <chowdsp_plugin_utils/chowdsp_plugin_utils.h>

// values oh ButterworthFilter 2 freq 7750 q 1.95
// values ch ButterworthFilter 2 freq 11450 q 1.55
// values ch output gain 2 FirstOrderHPF 8000 -> FirstOrderLPF 11000 ( -1500 <- 9500 -> +1500)
// values oh output gain 1.3 FirstOrderHPF 9500 -> FirstOrderLPF 12500 ( -1500 <- 11000 -> +1500)
// reso Resonator rfb 82e3 r_g 425 c 3.79e-9 gain 0.325


// closed etMin{-13.287}, etMax{2.0}; 	"tEnvA": 0.5e-3, D = 0, S = 0, "tEnvR": 60e-3, trigger = 1e-3

// #include "bar.hpp"
// #include "values_reso.hpp"
#include "Resonator.hpp"

void parseResponseSimulationFile(juce::String simFilename, std::vector<float> *destXVec, std::vector<float> *destYVec)
{
    const std::regex txt_regex("([0-9.eE+-]+)	\\(([0-9.eE+-]+)dB,(.*)\\).*");
    const juce::File simFile = juce::File(simFilename);
    juce::StringArray linesArray;
    simFile.readLines(linesArray);
    // printf("Length linesArray %u\n", linesArray.size());
    for (size_t i = 0; i < linesArray.size(); i++) {
        if (i == 0)
            continue;
        std::smatch base_match;
        std::string curLine = linesArray[i].toStdString();
        if (std::regex_match(curLine, base_match, txt_regex))
        {
            // printf("%f\n", juce::String(base_match[1].str()).getFloatValue());
            // printf("%f\n", juce::String(base_match[2].str()).getFloatValue());
            destXVec->push_back(juce::String(base_match[1].str()).getFloatValue());
            destYVec->push_back(juce::String(base_match[2].str()).getFloatValue());
        }
        // printf("%s -> %f %f\n", curLine.c_str(), juce::String(base_match[1].str()).getFloatValue(), juce::String(base_match[2].str()).getFloatValue());
    }
}

void parseDCSimulationFile(juce::String simFilename, std::vector<float> *destXVec, std::vector<float> *destYVec)
{
    const std::regex txt_regex("([0-9.eE+-]+)	([0-9.eE+-]+).*");
    const juce::File simFile = juce::File(simFilename);
    juce::StringArray linesArray;
    simFile.readLines(linesArray);
    // printf("Length linesArray %u\n", linesArray.size());
    for (size_t i = 0; i < linesArray.size(); i++) {
        if (i == 0)
            continue;
        std::smatch base_match;
        std::string curLine = linesArray[i].toStdString();
        if (std::regex_match(curLine, base_match, txt_regex))
        {
            // printf("%f\n", juce::String(base_match[1].str()).getFloatValue());
            // printf("%f\n", juce::String(base_match[2].str()).getFloatValue());
            destXVec->push_back(juce::String(base_match[1].str()).getFloatValue());
            destYVec->push_back(juce::String(base_match[2].str()).getFloatValue());
        }
        // printf("%s -> %f %f\n", curLine.c_str(), juce::String(base_match[1].str()).getFloatValue(), juce::String(base_match[2].str()).getFloatValue());
    }
}

chowdsp::GenericTweaksFile<false> options;

// plot response
void plotResponse()
{
    auto freq = options.getProperty<float> ("freq");
    auto freq2 = options.getProperty<float> ("freq2");
    auto freq3 = options.getProperty<float> ("freq3");
    auto q = options.getProperty<float> ("q");
    auto gain = options.getProperty<float> ("gain");
    auto xmin = options.getProperty<float> ("xmin");
    auto xmax = options.getProperty<float> ("xmax");
    auto ymin = options.getProperty<float> ("ymin");
    auto ymax = options.getProperty<float> ("ymax");
    auto Rfb = options.getProperty<float> ("Rfb");
    auto R_g = options.getProperty<float> ("R_g");
    auto C = options.getProperty<float> ("C");
    auto simFilename = options.getProperty<juce::String> ("file");

    std::vector<float> freqValues;
    std::vector<float> dbValues;
    parseResponseSimulationFile(simFilename, &freqValues, &dbValues);

    HatResonatorWDF reso;
    reso.prepare(96000);
    reso.setParameters(Rfb, R_g, C);
    chowdsp::Gain<float> preGain;
    preGain.setGainLinear(gain);

    chowdsp::FirstOrderHPF<float> fi;
    chowdsp::FirstOrderLPF<float> fi2;

    fi.reset();
    fi.calcCoefs(freq3 - 1500, 96000.0);
    fi2.reset();
    fi2.calcCoefs(freq3 + 1500, 96000.0);

    chowdsp::SpectrumPlotBase base {
            chowdsp::SpectrumPlotParams {
                500.0f,
                20000.0f,
                -20.0f,
                30.0f }
    };
    chowdsp::GenericFilterPlotter plotter { base, {96000.0f, 1.0f / 12.0f, 13} };
    plotter.runFilterCallback = [&fi, &fi2, &preGain, gain, freq, freq2, q, &reso] (const float* in, float* out, int N)
    {   
        std::copy (in, in + N, out);
        preGain.prepare({96000.0, juce::uint32(N), 1});
        preGain.process (chowdsp::BufferView { out, N });
        // reso.reset();
        // reso.prepare(96000);
        // reso.setParameters(8000.0, 0.06);
        // for (size_t i = 0; i < N; i++) {
        //     out[i] = reso.processSample(out[i]);
        // }
        fi.prepare({96000.0, juce::uint32(N), 1});
        fi.processBlock (chowdsp::BufferView { out, N });
        fi2.prepare({96000.0, juce::uint32(N), 1});
        fi2.processBlock (chowdsp::BufferView { out, N });
    };

    auto [freqAxis, magAxis] = plotter.plotFilterMagnitudeResponse();

    Plot2D plot;
    plot.xlabel("Frequency (Hz)");
    plot.ylabel("Magnitude (dB)");
    plot.drawCurve(freqAxis, magAxis).label("");
    plot.drawCurve(freqValues, dbValues).label("LTSpice");
    plot.xtics().logscale();
    plot.xrange(xmin, xmax);
    plot.yrange(ymin, ymax);
    Figure fig = {{plot}};
    Canvas canvas = {{fig}};
    canvas.size(1300,1300);
    canvas.show();
    // canvas.save("/tmp/plot.png")
}

struct EnvTime
{
    static constexpr float defaultEtMin{-8}, defaultEtMax{3.32192809489}; // log2(10)
    float etMin{defaultEtMin}, etMax{defaultEtMax};
    
    EnvTime() {}
    EnvTime(float etMin, float etMax)
        : etMin(etMin), etMax(etMax) {}
    
    float value{etMin};

    float getValue() { return value; }
    void setValue(float nvalue) { value = nvalue; }
    
    std::string getDisplayValueString()
    {
        auto v = getValue() * (etMax - etMin) + etMin;

        if (getValue() < 0.0001)
        {
            std::string mv;
            if (getMinString(mv))
            {
                return mv;
            }
        }
        char valstring[32];
        sprintf(valstring, "%.4f s", pow(2, v));
        return std::string(valstring);
    }
    
    void printDisplayValueString(float value)
    {
        auto v = value * (etMax - etMin) + etMin;
        char valstring[32];
        sprintf(valstring, "%.4f s", pow(2, v));
        printf("Value %f -> Time %s\n", value, valstring);

    }
    void setDisplayValueString(std::string s)
    {
        auto q = std::atof(s.c_str());
        auto v = log2(std::clamp(q, pow(2., etMin), pow(2., etMax)));
        auto vn = (v - etMin) / (etMax - etMin);
        setValue(vn);
    }
    
    float getValueTime(float q)
    {
        auto v = log2(std::clamp(double(q), pow(2., etMin), pow(2., etMax)));
        auto vn = (v - etMin) / (etMax - etMin);
        return vn;
    }

    virtual bool getMinString(std::string &s) { return false; }
};

static constexpr int tbs{8};
struct SRProvider
{
    double samplerate{96000}, sampleRateInv{1.f / samplerate};
    float envelope_rate_linear_nowrap(float f) const { return tbs * sampleRateInv * pow(2.f, -f); }
} srp;

struct ShortRange
{
    // 100µs -> 4s
    static constexpr float etMin{-13.287}, etMax{2.0};
};

void plotEnvelope()
{   
    auto xmin = options.getProperty<float> ("xmin");
    auto xmax = options.getProperty<float> ("xmax");
    auto ymin = options.getProperty<float> ("ymin");
    auto ymax = options.getProperty<float> ("ymax");
    auto simFilename = options.getProperty<juce::String> ("file");
    auto tEnvA = options.getProperty<float> ("tEnvA");
    auto tEnvR = options.getProperty<float> ("tEnvR");
    auto mult = options.getProperty<float> ("mult");
    auto trigDuration = options.getProperty<float> ("trigDuration");
    
    rack::dsp::PulseGenerator trigger;
    
    EnvTime et1(ShortRange::etMin, ShortRange::etMax);
    
    // et1.setDisplayValueString(std::string("0.1 s"));
    // printf("%s\n", et1.getDisplayValueString().c_str());
    
    et1.printDisplayValueString(0.0);
    et1.printDisplayValueString(1.0);
    
    auto adsr = sst::basic_blocks::modulators::ADSREnvelope<SRProvider, tbs, ShortRange>(&srp);
    
    float deltaTime = 1./96000.;
    float length = 0.2;
    
    std::vector<float> adsr_vec;
    std::vector<float> time_vec;
    // std::vector<float> out2;
    // std::vector<float> out3;
    // std::vector<float> time_x;
    
    float EnvA = et1.getValueTime(tEnvA);
    float EnvD = et1.getValueTime(0.2);
    float EnvS = 1.0;
    float EnvR = et1.getValueTime(tEnvR);
    
    
    float time = 0.0f;
    while (time < length) {
        if (time == int(0.001 / deltaTime) * deltaTime)
        {
            trigger.trigger(trigDuration);
            adsr.attackFrom(0.0, EnvA, 1, false); // initial, attacktime, ashp, digital?
        }
        auto b = trigger.process(deltaTime);
        adsr.process(EnvA, EnvD, EnvS, EnvR, 1, 1, 1, b); // a, d, s, r, ashp, dshp, rshp, gateActive
        adsr_vec.push_back(adsr.output * mult);
        time_vec.push_back(time);
        
        time += deltaTime;
    }
    
    std::vector<float> timeValues;
    std::vector<float> voltageValues;
    parseDCSimulationFile(simFilename, &timeValues, &voltageValues);
    
    Plot2D plot;
    plot.xlabel("Time (s)");
    plot.ylabel("Value");
    plot.drawCurve(time_vec, adsr_vec).label("ADSR");
    plot.drawCurve(timeValues, voltageValues).label("LTSpice");
    plot.xrange(xmin, xmax);
    plot.yrange(ymin, ymax);
    Figure fig = {{plot}};
    Canvas canvas = {{fig}};
    canvas.size(1300,1300);
    canvas.show();
    // canvas.save("/tmp/plot.png");
}


auto main(int argc, const char *argv[]) -> int
{   
    const juce::File optionsFileJ = juce::File(argv[1]);
    options.initialise (optionsFileJ, 1);
    
    // plotResponse();
    plotEnvelope();
    
    return 0;
}
