/*
 * vad.cc
 *
 *  Created on: Apr 2, 2018
 *      Author: tao
 */

#include <assert.h>
#include "matrix/Vector.h"
#include "ivector/vad.h"

void ComputeVadEnergy(const VadEnergyOptions &opts,
                      const Matrix &feats,
                      Vector *output_voiced) {
  int32 T = feats.NumRows();
  output_voiced->Resize(T);
  if (T == 0) {
    std::cout << "Empty features";
    return;
  }
  Vector log_energy(T);
  log_energy.CopyColFromMat(feats, 0); // column zero is log-energy.

  BaseFloat energy_threshold = opts.vad_energy_threshold;
  if (opts.vad_energy_mean_scale != 0.0) {
    assert(opts.vad_energy_mean_scale > 0.0);
    energy_threshold += opts.vad_energy_mean_scale * log_energy.Sum() / T;
  }

  assert(opts.vad_frames_context >= 0);
  assert(opts.vad_proportion_threshold > 0.0 &&
               opts.vad_proportion_threshold < 1.0);
  for (int32 t = 0; t < T; t++) {
    const BaseFloat *log_energy_data = log_energy.Data();
    int32 num_count = 0, den_count = 0, context = opts.vad_frames_context;
    for (int32 t2 = t - context; t2 <= t + context; t2++) {
      if (t2 >= 0 && t2 < T) {
        den_count++;
        if (log_energy_data[t2] > energy_threshold){
          num_count++;
        }
      }
    }
    if (num_count >= den_count * opts.vad_proportion_threshold)
      (*output_voiced)(t) = 1.0;
    else
      (*output_voiced)(t) = 0.0;
  }
}
