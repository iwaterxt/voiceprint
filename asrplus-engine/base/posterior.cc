/*
 * posterior.cc
 *
 *  Created on: Feb 18, 2018
 *      Author: tao
 */

#include "base/posterior.h"


BaseFloat VectorToPosteriorEntry(
  const Vector &log_likes,
  int32 num_gselect,
  BaseFloat min_post,
  std::vector<std::pair<int32, BaseFloat> > *post_entry) {
  assert(num_gselect > 0 && min_post >= 0 && min_post < 1.0);
  // we name num_gauss assuming each entry in log_likes represents a Gaussian;
  // it doesn't matter if they don't.
  int32 num_gauss = log_likes.Dim();
  assert(num_gauss > 0);
  if (num_gselect > num_gauss)
    num_gselect = num_gauss;
  Vector log_likes_normalized(log_likes);
  BaseFloat ans = log_likes_normalized.ApplySoftMax();
  std::vector<std::pair<int32, BaseFloat> > temp_post(num_gauss);
  for (int32 g = 0; g < num_gauss; g++)
    temp_post[g] = std::pair<int32, BaseFloat>(g, log_likes_normalized(g));
  CompareReverseSecond compare;
  // Sort in decreasing order on posterior.  For efficiency we
  // first do nth_element and then sort, as we only need the part we're
  // going to output, to be sorted.
  std::nth_element(temp_post.begin(),
                   temp_post.begin() + num_gselect, temp_post.end(),
                   compare);
  std::sort(temp_post.begin(), temp_post.begin() + num_gselect,
            compare);

  post_entry->clear();
  post_entry->insert(post_entry->end(),
                     temp_post.begin(), temp_post.begin() + num_gselect);
  while (post_entry->size() > 1 && post_entry->back().second < min_post)
    post_entry->pop_back();
  // Now renormalize to sum to one after pruning.
  BaseFloat tot = 0.0;
  size_t size = post_entry->size();
  for (size_t i = 0; i < size; i++)
    tot += (*post_entry)[i].second;
  BaseFloat inv_tot = 1.0 / tot;
  for (size_t i = 0; i < size; i++)
    (*post_entry)[i].second *= inv_tot;
  return ans;
}

void ScalePosterior(BaseFloat scale, Posterior *post) {
  if (scale == 1.0) return;
  for (size_t i = 0; i < post->size(); i++) {
    if (scale == 0.0) {
      (*post)[i].clear();
    } else {
      for (size_t j = 0; j < (*post)[i].size(); j++)
        (*post)[i][j].second *= scale;
    }
  }
}

BaseFloat TotalPosterior(const Posterior &post) {
  BaseFloat sum =  0.0;
  size_t T = post.size();
  for (size_t t = 0; t < T; t++) {
    size_t U = post[t].size();
    for (size_t i = 0; i < U; i++) {
      sum += post[t][i].second;
    }
  }
  return sum;
}
