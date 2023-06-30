#include <iostream>
#include <string>
#include <cstdio>

#include "lib.hpp"
#include "sst/basic-blocks/modulators/ADAREnvelope.h"

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

#include "digital.hpp"

static constexpr int tbs{16};
struct SampleSRProvider
{
    double samplerate{44100}, sampleRateInv{1.f / samplerate};
    float envelope_rate_linear_nowrap(float f) const { return tbs * sampleRateInv * pow(2.f, -f); }
} srp;

float out[44100*2];

rack::dsp::PulseGenerator trigger;

auto main() -> int
{
    auto const lib = library {};
    auto const message = "Hello from " + lib.name + "!";
    std::cout << message << '\n';
    
    drwav_data_format format;
    format.container = drwav_container_riff;     // <-- drwav_container_riff = normal WAV files, drwav_container_w64 = Sony Wave64.
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;          // <-- Any of the DR_WAVE_FORMAT_* codes.
    format.channels = 1;
    format.sampleRate = 44100;
    format.bitsPerSample = 32;
    drwav wav;
    drwav_init_file_write(&wav, "out.wav", &format, NULL);
    
    auto adar = sst::basic_blocks::modulators::ADAREnvelope<SampleSRProvider, tbs>(&srp);
    auto A = log2(1.5e-3f);
    auto H = log2(2.2e-3f);
    auto R = log2(100e-3f);
    
    
    float deltaTime = 1.0/44100.0;
    
    for (size_t s = 0; s < 44100*2; s++) {
        if (s == 200)
        {
            trigger.trigger(3.7e-3f);
            adar.attackFrom(0.0, 1, false, true); // initial, ashp, digital?, gated?
            // printf("Hey!\n");
        }
        auto b = trigger.process(deltaTime);
            
        adar.processScaledAD(A, R, 1, 1, b); // a, d, ashp, dshp, gateActive
        out[s] = adar.output * 2 - 1.0;
        // out[s] = b ? 1.0 : -1.0;
    }
    drwav_uint64 framesWritten = drwav_write_pcm_frames(&wav, 44100*2, out);
    drwav_uninit(&wav);
    
    return 0;
}
