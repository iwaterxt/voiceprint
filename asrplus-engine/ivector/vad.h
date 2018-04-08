/*
 * vad.h
 *
 *  Created on: Apr 2, 2018
 *      Author: tao
 */

#ifndef VAD_H_
#define VAD_H_

struct VadEnergyOptions {
  BaseFloat vad_energy_threshold;
  BaseFloat vad_energy_mean_scale;
  int32 vad_frames_context;
  BaseFloat vad_proportion_threshold;

  VadEnergyOptions(): vad_energy_threshold(5.5),
                      vad_energy_mean_scale(0.5),
                      vad_frames_context(0),
                      vad_proportion_threshold(0.6) { }
  void Register(OptionsItf *opts) {
    opts->Register("vad-energy-threshold", &vad_energy_threshold,
                   "Constant term in energy threshold for MFCC0 for VAD (also see "
                   "--vad-energy-mean-scale)");
    opts->Register("vad-energy-mean-scale", &vad_energy_mean_scale,
                   "If this is set to s, to get the actual threshold we "
                   "let m be the mean log-energy of the file, and use "
                   "s*m + vad-energy-threshold");
    opts->Register("vad-frames-context", &vad_frames_context,
                   "Number of frames of context on each side of central frame, "
                   "in window for which energy is monitored");
    opts->Register("vad-proportion-threshold", &vad_proportion_threshold,
                   "Parameter controlling the proportion of frames within "
                   "the window that need to have more energy than the "
                   "threshold");
  }
};

void ComputeVadEnergy(const VadEnergyOptions &opts,
                      const Matrix &input_features,
                      Vector *output_voiced);



#endif /* VAD_H_ */
