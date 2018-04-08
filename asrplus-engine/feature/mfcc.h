/*
 * mfcc.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#ifndef MFCC_H_
#define MFCC_H_

#include "feature/feature-common.h"
#include "matrix/srfft.h"
#include <map>

struct MfccOptions {
  FrameExtractionOptions frame_opts;
  MelBanksOptions mel_opts;
  int32 num_ceps;  // e.g. 13: num cepstral coeffs, counting zero.
  bool use_energy;  // use energy; else C0
  BaseFloat energy_floor;
  bool raw_energy;  // If true, compute energy before preemphasis and windowing
  BaseFloat cepstral_lifter;  // Scaling factor on cepstra for HTK compatibility.
                              // if 0.0, no liftering is done.
  bool htk_compat;  // if true, put energy/C0 last and introduce a factor of
                    // sqrt(2) on C0 to be the same as HTK.

  MfccOptions() : mel_opts(23),
                  // defaults the #mel-banks to 23 for the MFCC computations.
                  // this seems to be common for 16khz-sampled data,
                  // but for 8khz-sampled data, 15 may be better.
                  num_ceps(20),
                  use_energy(true),
                  energy_floor(0.0),  // not in log scale: a small value e.g. 1.0e-10
                  raw_energy(true),
                  cepstral_lifter(22.0),
                  htk_compat(false) {}

  void Register(OptionsItf *opts) {
    frame_opts.Register(opts);
    mel_opts.Register(opts);
    opts->Register("num-ceps", &num_ceps,
                   "Number of cepstra in MFCC computation (including C0)");
    opts->Register("use-energy", &use_energy,
                   "Use energy (not C0) in MFCC computation");
    opts->Register("energy-floor", &energy_floor,
                   "Floor on energy (absolute, not relative) in MFCC computation");
    opts->Register("raw-energy", &raw_energy,
                   "If true, compute energy before preemphasis and windowing");
    opts->Register("cepstral-lifter", &cepstral_lifter,
                   "Constant that controls scaling of MFCCs");
    opts->Register("htk-compat", &htk_compat,
                   "If true, put energy or C0 last and use a factor of sqrt(2) on "
                   "C0.  Warning: not sufficient to get HTK compatible features "
                   "(need to change other parameters).");
  }
};



// This is the new-style interface to the MFCC computation.
class MfccComputer {
 public:
  typedef MfccOptions Options;
  explicit MfccComputer(const MfccOptions &opts);
  MfccComputer(const MfccComputer &other);

  const FrameExtractionOptions &GetFrameOptions() const {
    return opts_.frame_opts;
  }

  int32 Dim() const { return opts_.num_ceps; }

  bool NeedRawLogEnergy() { return opts_.use_energy && opts_.raw_energy; }

  void Compute(BaseFloat signal_log_energy,
               BaseFloat vtln_warp,
               Vector *signal_frame,
               Vector *feature);

  ~MfccComputer();
 private:
  // disallow assignment.
  MfccComputer &operator = (const MfccComputer &in);

  const MelBanks *GetMelBanks(BaseFloat vtln_warp);

  MfccOptions opts_;
  Vector lifter_coeffs_;
  Matrix dct_matrix_;  // matrix we left-multiply by to perform DCT.
  BaseFloat log_energy_floor_;
  std::map<BaseFloat, MelBanks*> mel_banks_;  // BaseFloat is VTLN coefficient.
  SplitRadixRealFft *srfft_;

  // note: mel_energies_ is specific to the frame we're processing, it's
  // just a temporary workspace.
  Vector mel_energies_;
};

typedef OfflineFeatureTpl<MfccComputer> Mfcc;


#endif /* MFCC_H_ */
