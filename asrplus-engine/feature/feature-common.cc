/*
 * feature-common.cc
 *
 *  Created on: Feb 8, 2018
 *      Author: tao
 */

#include <assert.h>
#include <stdio.h>
#include "base/common.h"
#include "util/math.h"
#include "matrix/Vector.h"
#include "matrix/Matrix.h"
#include "feature/feature-common.h"
#include "feature/mfcc.h"
#include "feature/feature-function.h"

FeatureWindowFunction::FeatureWindowFunction(const FrameExtractionOptions &opts) {
  int32 frame_length = opts.WindowSize();
  assert(frame_length > 0);
  window.Resize(frame_length);
  double a = M_2PI / (frame_length-1);
  for (int32 i = 0; i < frame_length; i++) {
    double i_fl = static_cast<double>(i);
    if (opts.window_type == "hanning") {
      window(i) = 0.5  - 0.5*cos(a * i_fl);
    } else if (opts.window_type == "hamming") {
      window(i) = 0.54 - 0.46*cos(a * i_fl);
    } else if (opts.window_type == "povey") {  // like hamming but goes to zero at edges.
      window(i) = pow(0.5 - 0.5*cos(a * i_fl), 0.85);
    } else if (opts.window_type == "rectangular") {
      window(i) = 1.0;
    } else if (opts.window_type == "blackman") {
      window(i) = opts.blackman_coeff - 0.5*cos(a * i_fl) +
        (0.5 - opts.blackman_coeff) * cos(2 * a * i_fl);
    } else {
      std::cout << "Invalid window type " << opts.window_type<<std::endl;
    }
  }
}

void MelBanks::Compute(const Vector &power_spectrum,
                       Vector *mel_energies_out) const {
  int32 num_bins = bins_.size();
  assert(mel_energies_out->Dim() == num_bins);

  for (int32 i = 0; i < num_bins; i++) {
    int32 offset = bins_[i].first;
    const Vector &v(bins_[i].second);
    BaseFloat energy = VecVec(v, power_spectrum.Range(offset, v.Dim()));
    // HTK-like flooring- for testing purposes (we prefer dither)
    if (htk_mode_ && energy < 1.0) energy = 1.0;
    (*mel_energies_out)(i) = energy;

    // The following assert was added due to a problem with OpenBlas that
    // we had at one point (it was a bug in that library).  Just to detect
    // it early.
    //assert(!ISNAN((*mel_energies_out)(i)));
  }

}


MelBanks::MelBanks(const MelBanks &other):
    center_freqs_(other.center_freqs_),
    bins_(other.bins_),
    debug_(other.debug_),
    htk_mode_(other.htk_mode_) { }

MelBanks::MelBanks(const MelBanksOptions &opts,
                   const FrameExtractionOptions &frame_opts,
                   BaseFloat vtln_warp_factor):
    htk_mode_(opts.htk_mode) {
  int32 num_bins = opts.num_bins;
  if (num_bins < 3) std::cerr << "Must have at least 3 mel bins"<<std::endl;
  BaseFloat sample_freq = frame_opts.samp_freq;
  int32 window_length = static_cast<int32>(frame_opts.samp_freq*0.001*frame_opts.frame_length_ms);
  int32 window_length_padded =
      (frame_opts.round_to_power_of_two ?
       RoundUpToNearestPowerOfTwo(window_length) :
       window_length);
  assert(window_length_padded % 2 == 0);
  int32 num_fft_bins = window_length_padded/2;
  BaseFloat nyquist = 0.5 * sample_freq;

  BaseFloat low_freq = opts.low_freq, high_freq;
  if (opts.high_freq > 0.0)
    high_freq = opts.high_freq;
  else
    high_freq = nyquist + opts.high_freq;

  if (low_freq < 0.0 || low_freq >= nyquist
      || high_freq <= 0.0 || high_freq > nyquist
      || high_freq <= low_freq)
    std::cerr << "Bad values in options: low-freq " << low_freq
              << " and high-freq " << high_freq << " vs. nyquist "
              << nyquist<<std::endl;

  BaseFloat fft_bin_width = sample_freq / window_length_padded;
  // fft-bin width [think of it as Nyquist-freq / half-window-length]

  BaseFloat mel_low_freq = MelScale(low_freq);
  BaseFloat mel_high_freq = MelScale(high_freq);

  debug_ = opts.debug_mel;

  // divide by num_bins+1 in next line because of end-effects where the bins
  // spread out to the sides.
  BaseFloat mel_freq_delta = (mel_high_freq - mel_low_freq) / (num_bins+1);

  BaseFloat vtln_low = opts.vtln_low,
      vtln_high = opts.vtln_high;
  if (vtln_high < 0.0) vtln_high += nyquist;

  if (vtln_warp_factor != 1.0 &&
      (vtln_low < 0.0 || vtln_low <= low_freq
       || vtln_low >= high_freq
       || vtln_high <= 0.0 || vtln_high >= high_freq
       || vtln_high <= vtln_low))
    std::cerr << "Bad values in options: vtln-low " << vtln_low
              << " and vtln-high " << vtln_high << ", versus "
              << "low-freq " << low_freq << " and high-freq "
              << high_freq<<std::endl;

  bins_.resize(num_bins);
  center_freqs_.Resize(num_bins);

  for (int32 bin = 0; bin < num_bins; bin++) {
    BaseFloat left_mel = mel_low_freq + bin * mel_freq_delta,
        center_mel = mel_low_freq + (bin + 1) * mel_freq_delta,
        right_mel = mel_low_freq + (bin + 2) * mel_freq_delta;

    if (vtln_warp_factor != 1.0) {
      left_mel = VtlnWarpMelFreq(vtln_low, vtln_high, low_freq, high_freq,
                                 vtln_warp_factor, left_mel);
      center_mel = VtlnWarpMelFreq(vtln_low, vtln_high, low_freq, high_freq,
                                 vtln_warp_factor, center_mel);
      right_mel = VtlnWarpMelFreq(vtln_low, vtln_high, low_freq, high_freq,
                                  vtln_warp_factor, right_mel);
    }
    center_freqs_(bin) = InverseMelScale(center_mel);
    // this_bin will be a vector of coefficients that is only
    // nonzero where this mel bin is active.
    Vector this_bin(num_fft_bins);
    int32 first_index = -1, last_index = -1;
    for (int32 i = 0; i < num_fft_bins; i++) {
      BaseFloat freq = (fft_bin_width * i);  // Center frequency of this fft
                                             // bin.
      BaseFloat mel = MelScale(freq);
      if (mel > left_mel && mel < right_mel) {
        BaseFloat weight;
        if (mel <= center_mel)
          weight = (mel - left_mel) / (center_mel - left_mel);
        else
         weight = (right_mel-mel) / (right_mel-center_mel);
        this_bin(i) = weight;
        if (first_index == -1)
          first_index = i;
        last_index = i;
      }
    }
    assert(first_index != -1 && last_index >= first_index
                 && "You may have set --num-mel-bins too large.");

    bins_[bin].first = first_index;
    int32 size = last_index + 1 - first_index;
    bins_[bin].second.Resize(size);
    bins_[bin].second.CopyFromVec(this_bin.Range(first_index, size));

    // Replicate a bug in HTK, for testing purposes.
    if (opts.htk_mode && bin == 0 && mel_low_freq != 0.0)
      bins_[bin].second(0) = 0.0;

  }
}

BaseFloat MelBanks::VtlnWarpFreq(BaseFloat vtln_low_cutoff,  // upper+lower frequency cutoffs for VTLN.
                                 BaseFloat vtln_high_cutoff,
                                 BaseFloat low_freq,  // upper+lower frequency cutoffs in mel computation
                                 BaseFloat high_freq,
                                 BaseFloat vtln_warp_factor,
                                 BaseFloat freq) {

  if (freq < low_freq || freq > high_freq) return freq;  // in case this gets called
  // for out-of-range frequencies, just return the freq.

  assert(vtln_low_cutoff > low_freq &&
               "be sure to set the --vtln-low option higher than --low-freq");
  assert(vtln_high_cutoff < high_freq &&
               "be sure to set the --vtln-high option lower than --high-freq [or negative]");
  BaseFloat one = 1.0;
  BaseFloat l = vtln_low_cutoff * std::max(one, vtln_warp_factor);
  BaseFloat h = vtln_high_cutoff * std::min(one, vtln_warp_factor);
  BaseFloat scale = 1.0 / vtln_warp_factor;
  BaseFloat Fl = scale * l;  // F(l);
  BaseFloat Fh = scale * h;  // F(h);
  assert(l > low_freq && h < high_freq);
  // slope of left part of the 3-piece linear function
  BaseFloat scale_left = (Fl - low_freq) / (l - low_freq);
  // [slope of center part is just "scale"]

  // slope of right part of the 3-piece linear function
  BaseFloat scale_right = (high_freq - Fh) / (high_freq - h);

  if (freq < l) {
    return low_freq + scale_left * (freq - low_freq);
  } else if (freq < h) {
    return scale * freq;
  } else {  // freq >= h
    return high_freq + scale_right * (freq - high_freq);
  }
}


BaseFloat MelBanks::VtlnWarpMelFreq(BaseFloat vtln_low_cutoff,  // upper+lower frequency cutoffs for VTLN.
                                    BaseFloat vtln_high_cutoff,
                                    BaseFloat low_freq,  // upper+lower frequency cutoffs in mel computation
                                    BaseFloat high_freq,
                                    BaseFloat vtln_warp_factor,
                                    BaseFloat mel_freq) {
  return MelScale(VtlnWarpFreq(vtln_low_cutoff, vtln_high_cutoff,
                               low_freq, high_freq,
                               vtln_warp_factor, InverseMelScale(mel_freq)));
}



template <class F>
void OfflineFeatureTpl<F>::Compute(
    const Vector &wave,
    BaseFloat vtln_warp,
    Matrix *output) {
  assert(output != NULL);
  int32 rows_out = NumFrames(wave.Dim(), computer_.GetFrameOptions()),
      cols_out = computer_.Dim();
  if (rows_out == 0) {
    output->Resize(0, 0);
    return;
  }
  output->Resize(rows_out, cols_out);
  Vector window;  // windowed waveform.
  bool use_raw_log_energy = computer_.NeedRawLogEnergy();
  for (int32 r = 0; r < rows_out; r++) {  // r is frame index.
    BaseFloat raw_log_energy = 0.0;
    ExtractWindow(0, wave, r, computer_.GetFrameOptions(),
                  feature_window_function_, &window,
                  (use_raw_log_energy ? &raw_log_energy : NULL));
    Vector output_row(*output, r);
    computer_.Compute(raw_log_energy, vtln_warp, &window, &output_row);
    output->CopyRowFromVec(output_row, r);
  }
}

template <class F>
void OfflineFeatureTpl<F>::Compute(
    const Vector &wave,
    BaseFloat vtln_warp,
    Matrix *output) const {
  OfflineFeatureTpl<F> temp(*this);
  // call the non-const version of Compute() on a temporary copy of this object.
  // This is a workaround for const-ness that may sometimes be useful in
  // multi-threaded code, although it's not optimally efficient.
  temp.Compute(wave, vtln_warp, output);
}

template class OfflineFeatureTpl<MfccComputer>;
