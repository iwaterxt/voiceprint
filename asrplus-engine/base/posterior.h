/*
 * posterior.h
 *
 *  Created on: Oct 10, 2017
 *      Author: tao
 */

#ifndef POSTERIOR_H_
#define POSTERIOR_H_

#include <vector>
#include <algorithm>
#include <assert.h>
#include <utility>
#include "matrix/Vector.h"

typedef std::vector<std::vector<std::pair<int32, BaseFloat> > > Posterior;

struct CompareReverseSecond {
  // view this as an "<" operator used for sorting, except it behaves like
  // a ">" operator on the .second field of the pair because we want the
  // sort to be in reverse order (greatest to least) on posterior.
  bool operator() (const std::pair<int32, BaseFloat> &a,
                   const std::pair<int32, BaseFloat> &b) {
    return (a.second > b.second);
  }
};


BaseFloat VectorToPosteriorEntry(
  const Vector &log_likes,
  int32 num_gselect,
  BaseFloat min_post,
  std::vector<std::pair<int32, BaseFloat> > *post_entry) ;

void ScalePosterior(BaseFloat scale, Posterior *post);

BaseFloat TotalPosterior(const Posterior &post);
#endif /* POSTERIOR_H_ */
