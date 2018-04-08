/*
 * feature-function.cc
 *
 *  Created on: Feb 7, 2018
 *      Author: tao
 */
#include "feature/feature-function.h"
#include "matrix/Vector.h"
#include "util/math.h"
#include <assert.h>

void ComputePowerSpectrum(Vector *waveform) {
  int32 dim = waveform->Dim();

  // no, letting it be non-power-of-two for now.
  // KALDI_ASSERT(dim > 0 && (dim & (dim-1) == 0));  // make sure a power of two.. actually my FFT code
  // does not require this (dan) but this is better in case we use different code [dan].

  // RealFft(waveform, true);  // true == forward (not inverse) FFT; makes no difference here,
  // as we just want power spectrum.

  // now we have in waveform, first half of complex spectrum
  // it's stored as [real0, realN/2-1, real1, im1, real2, im2, ...]
  int32 half_dim = dim/2;
  double first_energy = (*waveform)(0) * (*waveform)(0),
      last_energy = (*waveform)(1) * (*waveform)(1);  // handle this special case
  for (int32 i = 1; i < half_dim; i++) {
    double real = (*waveform)(i*2), im = (*waveform)(i*2 + 1);
    (*waveform)(i) = real*real + im*im;
  }
  (*waveform)(0) = first_energy;
  (*waveform)(half_dim) = last_energy;  // Will actually never be used, and anyway
  // if the signal has been bandlimited sensibly this should be zero.
}

void Dither(Vector *waveform, BaseFloat dither_value) {
  if (dither_value == 0.0)
    return;
  int32 dim = waveform->Dim();
  BaseFloat *data = waveform->Data();
  RandomState rstate;
  for (int32 i = 0; i < dim; i++)
    data[i] += RandGauss(&rstate) * dither_value;
}

int32 NumFrames(int64 num_samples,
                const FrameExtractionOptions &opts,
                bool flush) {
  int64 frame_shift = opts.WindowShift();
  int64 frame_length = opts.WindowSize();
  if (opts.snip_edges) {
    // with --snip-edges=true (the default), we use a HTK-like approach to
    // determining the number of frames-- all frames have to fit completely into
    // the waveform, and the first frame begins at sample zero.
    if (num_samples < frame_length)
      return 0;
    else
      return (1 + ((num_samples - frame_length) / frame_shift));
    // You can understand the expression above as follows: 'num_samples -
    // frame_length' is how much room we have to shift the frame within the
    // waveform; 'frame_shift' is how much we shift it each time; and the ratio
    // is how many times we can shift it (integer arithmetic rounds down).
  } else {
    // if --snip-edges=false, the number of frames is determined by rounding the
    // (file-length / frame-shift) to the nearest integer.  The point of this
    // formula is to make the number of frames an obvious and predictable
    // function of the frame shift and signal length, which makes many
    // segmentation-related questions simpler.
    //
    // Because integer division in C++ rounds toward zero, we add (half the
    // frame-shift minus epsilon) before dividing, to have the effect of
    // rounding towards the closest integer.
    int32 num_frames = (num_samples + (frame_shift / 2)) / frame_shift;

    if (flush)
      return num_frames;

    // note: 'end' always means the last plus one, i.e. one past the last.
    int64 end_sample_of_last_frame = FirstSampleOfFrame(num_frames - 1, opts)
        + frame_length;

    // the following code is optimized more for clarity than efficiency.
    // If flush == false, we can't output frames that extend past the end
    // of the signal.
    while (num_frames > 0 && end_sample_of_last_frame > num_samples) {
      num_frames--;
      end_sample_of_last_frame -= frame_shift;
    }
    return num_frames;
  }
}


int64 FirstSampleOfFrame(int32 frame,
                         const FrameExtractionOptions &opts) {
  int64 frame_shift = opts.WindowShift();
  if (opts.snip_edges) {
    return frame * frame_shift;
  } else {
    int64 midpoint_of_frame = frame_shift * frame  +  frame_shift / 2,
        beginning_of_frame = midpoint_of_frame  -  opts.WindowSize() / 2;
    return beginning_of_frame;
  }
}

void Preemphasize(Vector *waveform, BaseFloat preemph_coeff) {
  if (preemph_coeff == 0.0) return;
  assert(preemph_coeff >= 0.0 && preemph_coeff <= 1.0);
  for (int32 i = waveform->Dim()-1; i > 0; i--)
    (*waveform)(i) -= preemph_coeff * (*waveform)(i-1);
  (*waveform)(0) -= preemph_coeff * (*waveform)(0);
}

void ProcessWindow(const FrameExtractionOptions &opts,
                   const FeatureWindowFunction &window_function,
                   Vector *window,
                   BaseFloat *log_energy_pre_window) {

  int32 frame_length = opts.WindowSize();
  assert(window->Dim() == frame_length);
  if (opts.dither != 0.0)
    Dither(window, opts.dither);
  if (opts.remove_dc_offset)
    window->Add(-window->Sum() / frame_length);
  if (log_energy_pre_window != NULL) {
    BaseFloat energy = std::max(VecVec(*window, *window),
                                std::numeric_limits<BaseFloat>::epsilon());
    *log_energy_pre_window = Log(energy);
  }
  if (opts.preemph_coeff != 0.0)
    Preemphasize(window, opts.preemph_coeff);
  window->MulElements(window_function.window);

}



void ExtractWindow(int64 sample_offset,
                   const Vector &wave,
                   int32 f,  // with 0 <= f < NumFrames(feats, opts)
                   const FrameExtractionOptions &opts,
                   const FeatureWindowFunction &window_function,
                   Vector *window,
                   BaseFloat *log_energy_pre_window) {
  assert(sample_offset >= 0 && wave.Dim() != 0);
  int32 frame_length = opts.WindowSize(),
      frame_length_padded = opts.PaddedWindowSize();
  int64 num_samples = sample_offset + wave.Dim(),
      start_sample = FirstSampleOfFrame(f, opts),
      end_sample = start_sample + frame_length;

  if (opts.snip_edges) {
    assert(start_sample >= sample_offset &&
                 end_sample <= num_samples);
  } else {
    assert(sample_offset == 0 || start_sample >= sample_offset);
  }

  if (window->Dim() != frame_length_padded)
    window->Resize(frame_length_padded);

  // wave_start and wave_end are start and end indexes into 'wave', for the
  // piece of wave that we're trying to extract.
  int32 wave_start = int32(start_sample - sample_offset),
      wave_end = wave_start + frame_length;
  if (wave_start >= 0 && wave_end <= wave.Dim()) {
    // the normal case-- no edge effects to consider.
    window->CopyFromVecRange(wave.Range(wave_start, frame_length));
  } else {
    // Deal with any end effects by reflection, if needed.  This code will only
    // be reached for about two frames per utterance, so we don't concern
    // ourselves excessively with efficiency.
    int32 wave_dim = wave.Dim();
    for (int32 s = 0; s < frame_length; s++) {
      int32 s_in_wave = s + wave_start;
      while (s_in_wave < 0 || s_in_wave >= wave_dim) {
        // reflect around the beginning or end of the wave.
        // e.g. -1 -> 0, -2 -> 1.
        // dim -> dim - 1, dim + 1 -> dim - 2.
        // the code supports repeated reflections, although this
        // would only be needed in pathological cases.
        if (s_in_wave < 0) s_in_wave = - s_in_wave - 1;
        else s_in_wave = 2 * wave_dim - 1 - s_in_wave;
      }
      (*window)(s) = wave(s_in_wave);
    }
  }

//  if (frame_length_padded > frame_length)
//    window->Range(frame_length, frame_length_padded - frame_length).SetZero();
  Vector frame(frame_length);
  frame.CopyFromData(window->Data(), frame_length);
  window->SetZero();

  ProcessWindow(opts, window_function, &frame, log_energy_pre_window);
  window->CopyFromVecRange(frame, frame_length);

}


