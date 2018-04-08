/*
 * mfcc.cc
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */
#include <limits>
#include <assert.h>
#include <map>
#include "base/common.h"
#include "matrix/srfft.h"
#include "matrix/Matrix.h"
#include "matrix/Vector.h"
#include "feature/mfcc.h"
#include "feature/feature-common.h"
#include "feature/feature-function.h"

void MfccComputer::Compute(BaseFloat signal_log_energy,
                           BaseFloat vtln_warp,
                           Vector *signal_frame,
                           Vector *feature) {
  assert(signal_frame->Dim() == opts_.frame_opts.PaddedWindowSize() &&
               feature->Dim() == this->Dim());

  const MelBanks &mel_banks = *(GetMelBanks(vtln_warp));

  if (opts_.use_energy && !opts_.raw_energy)
    signal_log_energy = Log(std::max(VecVec(*signal_frame, *signal_frame),
                                     std::numeric_limits<BaseFloat>::min()));
  if (srfft_ != NULL){  // Compute FFT using the split-radix algorithm.
    srfft_->Compute(signal_frame->Data(), true);
  }
  else  // An alternative algorithm that works for non-powers-of-two.
    RealFft(signal_frame, true);

  // Convert the FFT into a power spectrum.
  ComputePowerSpectrum(signal_frame);
  int32 dim = signal_frame->Dim() / 2 + 1 ;
  Vector power_spectrum(dim);
  power_spectrum.CopyFromData(signal_frame->Data(), signal_frame->Dim() / 2 + 1);

  mel_banks.Compute(power_spectrum, &mel_energies_);
  // avoid log of zero (which should be prevented anyway by dithering).
  mel_energies_.ApplyFloor(std::numeric_limits<BaseFloat>::epsilon());
  mel_energies_.ApplyLog();  // take the log.

  feature->SetZero();  // in case there were NaNs.
  // feature = dct_matrix_ * mel_energies [which now have log]
  feature->AddMatVec(1.0, dct_matrix_, kNoTrans, mel_energies_, 0.0);

  if (opts_.cepstral_lifter != 0.0)
    feature->MulElements(lifter_coeffs_);

  if (opts_.use_energy) {
    if (opts_.energy_floor > 0.0 && signal_log_energy < log_energy_floor_)
      signal_log_energy = log_energy_floor_;
    (*feature)(0) = signal_log_energy;
  }


}

MfccComputer::MfccComputer(const MfccOptions &opts):
    opts_(opts), srfft_(NULL),
    mel_energies_(opts.mel_opts.num_bins) {
  int32 num_bins = opts.mel_opts.num_bins;
  Matrix dct_matrix(num_bins, num_bins);
  ComputeDctMatrix(&dct_matrix);
  // Note that we include zeroth dct in either case.  If using the
  // energy we replace this with the energy.  This means a different
  // ordering of features than HTK.
  dct_matrix_.Resize(opts.num_ceps, num_bins);
  dct_matrix_.CopyFromMat(dct_matrix, 0, opts.num_ceps, 0, num_bins);
  if (opts.cepstral_lifter != 0.0) {
    lifter_coeffs_.Resize(opts.num_ceps);
    ComputeLifterCoeffs(opts.cepstral_lifter, &lifter_coeffs_);
  }
  if (opts.energy_floor > 0.0)
    log_energy_floor_ = Log(opts.energy_floor);

  int32 padded_window_size = opts.frame_opts.PaddedWindowSize();
  if ((padded_window_size & (padded_window_size-1)) == 0)  // Is a power of two...
    srfft_ = new SplitRadixRealFft(padded_window_size);

  // We'll definitely need the filterbanks info for VTLN warping factor 1.0.
  // [note: this call caches it.]
  GetMelBanks(1.0);
}

MfccComputer::MfccComputer(const MfccComputer &other):
    opts_(other.opts_), lifter_coeffs_(other.lifter_coeffs_),
    dct_matrix_(other.dct_matrix_),
    log_energy_floor_(other.log_energy_floor_),
    mel_banks_(other.mel_banks_),
    srfft_(NULL),
    mel_energies_(other.mel_energies_.Dim()) {
  for (std::map<BaseFloat, MelBanks*>::iterator iter = mel_banks_.begin();
       iter != mel_banks_.end(); ++iter)
    iter->second = new MelBanks(*(iter->second));
  if (other.srfft_ != NULL)
    srfft_ = new SplitRadixRealFft(*(other.srfft_));
}



MfccComputer::~MfccComputer() {
  for (std::map<BaseFloat, MelBanks*>::iterator iter = mel_banks_.begin();
      iter != mel_banks_.end();
      ++iter)
    delete iter->second;
  delete srfft_;
}

const MelBanks *MfccComputer::GetMelBanks(BaseFloat vtln_warp) {
  MelBanks *this_mel_banks = NULL;
  std::map<BaseFloat, MelBanks*>::iterator iter = mel_banks_.find(vtln_warp);
  if (iter == mel_banks_.end()) {
    this_mel_banks = new MelBanks(opts_.mel_opts,
                                  opts_.frame_opts,
                                  vtln_warp);
    mel_banks_[vtln_warp] = this_mel_banks;
  } else {
    this_mel_banks = iter->second;
  }
  return this_mel_banks;
}




