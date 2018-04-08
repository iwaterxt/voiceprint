/*
 * feature-function.h
 *
 *  Created on: Feb 7, 2018
 *      Author: tao
 */

#ifndef FEATURE_FUNCTION_H_
#define FEATURE_FUNCTION_H_
#include "matrix/Vector.h"
#include "feature/feature-common.h"

void ComputePowerSpectrum(Vector *complex_fft);

void Dither(Vector *waveform, BaseFloat dither_value);

int32 NumFrames(int64 num_samples,
                const FrameExtractionOptions &opts,
                bool flush = true);

int64 FirstSampleOfFrame(int32 frame,
                         const FrameExtractionOptions &opts);

void Preemphasize(Vector *waveform, BaseFloat preemph_coeff);

void ProcessWindow(const FrameExtractionOptions &opts,
                   const FeatureWindowFunction &window_function,
                   Vector *window,
                   BaseFloat *log_energy_pre_window = NULL);

void ExtractWindow(int64 sample_offset,
                   const Vector &wave,
                   int32 f,
                   const FrameExtractionOptions &opts,
                   const FeatureWindowFunction &window_function,
                   Vector *window,
                   BaseFloat *log_energy_pre_window = NULL);
#endif /* FEATURE_FUNCTION_H_ */
