/*
 * feature-common.h
 *
 *  Created on: Oct 10, 2017
 *      Author: tao
 */

#ifndef FEATURE_COMMON_H_
#define FEATURE_COMMON_H_

#include "util/math.h"


struct FrameExtractionOptions;

struct FeatureWindowFunction {
  FeatureWindowFunction() {}
  explicit FeatureWindowFunction(const FrameExtractionOptions &opts);
  FeatureWindowFunction(const FeatureWindowFunction &other):
      window(other.window) { }
  Vector window;
};



struct FrameExtractionOptions {
  BaseFloat samp_freq;
  BaseFloat frame_shift_ms;  // in milliseconds.
  BaseFloat frame_length_ms;  // in milliseconds.
  BaseFloat dither;  // Amount of dithering, 0.0 means no dither.
  BaseFloat preemph_coeff;  // Preemphasis coefficient.
  bool remove_dc_offset;  // Subtract mean of wave before FFT.
  std::string window_type;  // e.g. Hamming window
  bool round_to_power_of_two;
  BaseFloat blackman_coeff;
  bool snip_edges;
  bool allow_downsample;
  // May be "hamming", "rectangular", "povey", "hanning", "blackman"
  // "povey" is a window I made to be similar to Hamming but to go to zero at the
  // edges, it's pow((0.5 - 0.5*cos(n/N*2*pi)), 0.85)
  // I just don't think the Hamming window makes sense as a windowing function.
  FrameExtractionOptions():
      samp_freq(16000),
      frame_shift_ms(10.0),
      frame_length_ms(25.0),
      dither(0.0),
      preemph_coeff(0.97),
      remove_dc_offset(true),
      window_type("povey"),
      round_to_power_of_two(true),
      blackman_coeff(0.42),
      snip_edges(false),
      allow_downsample(false) { }

  void Register(OptionsItf *opts) {
    opts->Register("sample-frequency", &samp_freq,
                   "Waveform data sample frequency (must match the waveform file, "
                   "if specified there)");
    opts->Register("frame-length", &frame_length_ms, "Frame length in milliseconds");
    opts->Register("frame-shift", &frame_shift_ms, "Frame shift in milliseconds");
    opts->Register("preemphasis-coefficient", &preemph_coeff,
                   "Coefficient for use in signal preemphasis");
    opts->Register("remove-dc-offset", &remove_dc_offset,
                   "Subtract mean from waveform on each frame");
    opts->Register("dither", &dither, "Dithering constant (0.0 means no dither)");
    opts->Register("window-type", &window_type, "Type of window "
                   "(\"hamming\"|\"hanning\"|\"povey\"|\"rectangular\""
                   "|\"blackmann\")");
    opts->Register("blackman-coeff", &blackman_coeff,
                   "Constant coefficient for generalized Blackman window.");
    opts->Register("round-to-power-of-two", &round_to_power_of_two,
                   "If true, round window size to power of two by zero-padding "
                   "input to FFT.");
    opts->Register("snip-edges", &snip_edges,
                   "If true, end effects will be handled by outputting only frames that "
                   "completely fit in the file, and the number of frames depends on the "
                   "frame-length.  If false, the number of frames depends only on the "
                   "frame-shift, and we reflect the data at the ends.");
    opts->Register("allow-downsample", &allow_downsample,
                   "If true, allow the input waveform to have a higher frequency than"
                   "the specified --sample-frequency (and we'll downsample).");
  }
  int32 WindowShift() const {
    return static_cast<int32>(samp_freq * 0.001 * frame_shift_ms);
  }
  int32 WindowSize() const {
    return static_cast<int32>(samp_freq * 0.001 * frame_length_ms);
  }
  int32 PaddedWindowSize() const {
    return (round_to_power_of_two ? RoundUpToNearestPowerOfTwo(WindowSize()) :
                                    WindowSize());
  }
};

struct MelBanksOptions {
  int32 num_bins;  // e.g. 25; number of triangular bins
  BaseFloat low_freq;  // e.g. 20; lower frequency cutoff
  BaseFloat high_freq;  // an upper frequency cutoff; 0 -> no cutoff, negative
  // ->added to the Nyquist frequency to get the cutoff.
  BaseFloat vtln_low;  // vtln lower cutoff of warping function.
  BaseFloat vtln_high;  // vtln upper cutoff of warping function: if negative, added
                        // to the Nyquist frequency to get the cutoff.
  bool debug_mel;
  // htk_mode is a "hidden" config, it does not show up on command line.
  // Enables more exact compatibibility with HTK, for testing purposes.  Affects
  // mel-energy flooring and reproduces a bug in HTK.
  bool htk_mode;
  explicit MelBanksOptions(int num_bins = 25)
      : num_bins(num_bins), low_freq(20), high_freq(7800), vtln_low(100),
        vtln_high(-500), debug_mel(false), htk_mode(false) {}

  void Register(OptionsItf *opts) {
    opts->Register("num-mel-bins", &num_bins,
                   "Number of triangular mel-frequency bins");
    opts->Register("low-freq", &low_freq,
                   "Low cutoff frequency for mel bins");
    opts->Register("high-freq", &high_freq,
                   "High cutoff frequency for mel bins (if < 0, offset from Nyquist)");
    opts->Register("vtln-low", &vtln_low,
                   "Low inflection point in piecewise linear VTLN warping function");
    opts->Register("vtln-high", &vtln_high,
                   "High inflection point in piecewise linear VTLN warping function"
                   " (if negative, offset from high-mel-freq");
    opts->Register("debug-mel", &debug_mel,
                   "Print out debugging information for mel bin computation");
  }
};

class MelBanks {
 public:

  static inline BaseFloat InverseMelScale(BaseFloat mel_freq) {
    return 700.0f * (expf (mel_freq / 1127.0f) - 1.0f);
  }

  static inline BaseFloat MelScale(BaseFloat freq) {
    return 1127.0f * logf (1.0f + freq / 700.0f);
  }

  static BaseFloat VtlnWarpFreq(BaseFloat vtln_low_cutoff,
                                BaseFloat vtln_high_cutoff,  // discontinuities in warp func
                                BaseFloat low_freq,
                                BaseFloat high_freq,  // upper+lower frequency cutoffs in
                                // the mel computation
                                BaseFloat vtln_warp_factor,
                                BaseFloat freq);

  static BaseFloat VtlnWarpMelFreq(BaseFloat vtln_low_cutoff,
                                   BaseFloat vtln_high_cutoff,
                                   BaseFloat low_freq,
                                   BaseFloat high_freq,
                                   BaseFloat vtln_warp_factor,
                                   BaseFloat mel_freq);


  MelBanks(const MelBanksOptions &opts,
           const FrameExtractionOptions &frame_opts,
           BaseFloat vtln_warp_factor);

  /// Compute Mel energies (note: not log enerties).
  /// At input, "fft_energies" contains the FFT energies (not log).
  void Compute(const Vector &fft_energies,
               Vector *mel_energies_out) const;

  int32 NumBins() const { return bins_.size(); }

  // returns vector of central freq of each bin; needed by plp code.
  const Vector &GetCenterFreqs() const { return center_freqs_; }

  // Copy constructor
  MelBanks(const MelBanks &other);
 private:
  // Disallow assignment
  MelBanks &operator = (const MelBanks &other);

  // center frequencies of bins, numbered from 0 ... num_bins-1.
  // Needed by GetCenterFreqs().
  Vector center_freqs_;

  // the "bins_" vector is a vector, one for each bin, of a pair:
  // (the first nonzero fft-bin), (the vector of weights).
  std::vector<std::pair<int32, Vector > > bins_;

  bool debug_;
  bool htk_mode_;
};

template <class F>
class OfflineFeatureTpl {
 public:
  typedef typename F::Options Options;

  // Note: feature_window_function_ is the windowing function, which initialized
  // using the options class, that we cache at this level.
  OfflineFeatureTpl(const Options &opts):
      computer_(opts),
      feature_window_function_(computer_.GetFrameOptions()) { }

  // Internal (and back-compatibility) interface for computing features, which
  // requires that the user has already checked that the sampling frequency
  // of the waveform is equal to the sampling frequency specified in
  // the frame-extraction options.

  void Compute(const Vector &wave,
               BaseFloat vtln_warp,
               Matrix *output);

  // This const version of Compute() is a wrapper that
  // calls the non-const version on a temporary object.
  // It's less efficient than the non-const version.
  void Compute(const Vector &wave,
               BaseFloat vtln_warp,
               Matrix *output) const;


  int32 Dim() const { return computer_.Dim(); }

  // Copy constructor.
  OfflineFeatureTpl(const OfflineFeatureTpl<F> &other):
      computer_(other.computer_),
      feature_window_function_(other.feature_window_function_) { }
  private:
  // Disallow assignment.
  OfflineFeatureTpl<F> &operator =(const OfflineFeatureTpl<F> &other);

  F computer_;
  FeatureWindowFunction feature_window_function_;
};


#endif /* FEATURE_COMMON_H_ */
