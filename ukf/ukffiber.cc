/*
 * \file fiber.cc
 * \brief Implements functions defined in fiber.h
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#include "ukffiber.h"
#include <iostream>

void PostProcessFibers(const std::vector<UKFFiber> &raw_primary,
                       const std::vector<unsigned char> &discarded_fibers,
                       std::vector<UKFFiber> &fibers)
{
  assert(fibers.empty());
  const int num_half_fibers = static_cast<int>(raw_primary.size());

  assert((num_half_fibers > 0) && (num_half_fibers % 2 == 0));
  // if Noddi model is used Kappa is stored in trace, Vic in fa and Viso in freewater
  const bool record_rtop1 = !raw_primary[0].rtop1.empty();
  const bool record_rtop2 = !raw_primary[0].rtop2.empty();
  const bool record_rtop3 = !raw_primary[0].rtop3.empty();
  const bool record_rtop_model = !raw_primary[0].rtop_model.empty();
  const bool record_rtop_signal = !raw_primary[0].rtop_signal.empty();
  const bool record_free_water = !raw_primary[0].free_water.empty();
  const bool record_normMSE = !raw_primary[0].normMSE.empty();
  const bool record_cov = !raw_primary[0].covariance.empty();
  const bool record_w1 = !raw_primary[0].w1.empty();
  const bool record_w2 = !raw_primary[0].w2.empty();
  const bool record_w3 = !raw_primary[0].w3.empty();
  const bool record_w1w2angle = !raw_primary[0].w1w2angle.empty();
  const bool record_w1w3angle = !raw_primary[0].w1w3angle.empty();
  const bool record_Fm1 = !raw_primary[0].Fm1.empty();
  const bool record_lmd1 = !raw_primary[0].lmd1.empty();
  const bool record_Fm2 = !raw_primary[0].Fm2.empty();
  const bool record_lmd2 = !raw_primary[0].lmd2.empty();
  const bool record_Fm3 = !raw_primary[0].Fm3.empty();
  const bool record_lmd3 = !raw_primary[0].lmd3.empty();
  const bool record_varW1 = !raw_primary[0].varW1.empty();
  const bool record_varW2 = !raw_primary[0].varW2.empty();
  const bool record_varW3 = !raw_primary[0].varW3.empty();
  const bool record_varWiso = !raw_primary[0].varWiso.empty();
  const bool record_state = !raw_primary[0].state.empty();

  const int num_primary_fibers = num_half_fibers / 2;

  std::vector<int> num_points_on_primary_fiber(num_primary_fibers); // Number of points on each full primary fiber

  // Compute the numbers of points on full primary fibers
  for (int i = 0; i < num_primary_fibers; i++)
  {
    num_points_on_primary_fiber[i] =
        static_cast<int>(raw_primary[2 * i].position.size() + raw_primary[2 * i + 1].position.size()) - 1;
    // The two fibers share the same seed point
  }

  // Compute the number of valid full primary fibers and the number of valid full branches
  int num_valid_primary_fibers = 0;

  for (int i = 0; i < num_primary_fibers; i++)
  {
    if (num_points_on_primary_fiber[i] >= MINIMUM_NUM_POINTS_ON_FIBER && discarded_fibers[2 * i] != 1 && discarded_fibers[2 * i + 1] != 1)
      num_valid_primary_fibers++;
  }

  fibers.resize(num_valid_primary_fibers);

  // Join half fibers to full fibers
  int counter = 0;
  for (int i = 0; i < num_primary_fibers; i++)
  {
    if (num_points_on_primary_fiber[i] < MINIMUM_NUM_POINTS_ON_FIBER || discarded_fibers[2 * i] == 1 || discarded_fibers[2 * i + 1] == 1)
      continue;

    const UKFFiber &first_half = raw_primary[2 * i];
    const UKFFiber &second_half = raw_primary[2 * i + 1];

    fibers[counter].position.resize(num_points_on_primary_fiber[i]);
    if (record_rtop1)
      fibers[counter].rtop1.resize(num_points_on_primary_fiber[i]);
    if (record_rtop2)
      fibers[counter].rtop2.resize(num_points_on_primary_fiber[i]);
    if (record_rtop3)
      fibers[counter].rtop3.resize(num_points_on_primary_fiber[i]);
    if (record_rtop_model)
      fibers[counter].rtop_model.resize(num_points_on_primary_fiber[i]);
    if (record_rtop_signal)
      fibers[counter].rtop_signal.resize(num_points_on_primary_fiber[i]);
    if (record_free_water)
      fibers[counter].free_water.resize(num_points_on_primary_fiber[i]);
    if (record_w1)
      fibers[counter].w1.resize(num_points_on_primary_fiber[i]);
    if (record_w2)
      fibers[counter].w2.resize(num_points_on_primary_fiber[i]);
    if (record_w3)
      fibers[counter].w3.resize(num_points_on_primary_fiber[i]);
    if (record_w1w2angle)
      fibers[counter].w1w2angle.resize(num_points_on_primary_fiber[i]);
    if (record_w1w3angle)
      fibers[counter].w1w3angle.resize(num_points_on_primary_fiber[i]);
    if (record_Fm1)
      fibers[counter].Fm1.resize(num_points_on_primary_fiber[i]);
    if (record_lmd1)
      fibers[counter].lmd1.resize(num_points_on_primary_fiber[i]);
    if (record_Fm2)
      fibers[counter].Fm2.resize(num_points_on_primary_fiber[i]);
    if (record_lmd2)
      fibers[counter].lmd2.resize(num_points_on_primary_fiber[i]);
    if (record_Fm3)
      fibers[counter].Fm3.resize(num_points_on_primary_fiber[i]);
    if (record_lmd3)
      fibers[counter].lmd3.resize(num_points_on_primary_fiber[i]);
    if (record_varW1)
      fibers[counter].varW1.resize(num_points_on_primary_fiber[i]);
    if (record_varW2)
      fibers[counter].varW2.resize(num_points_on_primary_fiber[i]);
    if (record_varW3)
      fibers[counter].varW3.resize(num_points_on_primary_fiber[i]);
    if (record_varWiso)
      fibers[counter].varWiso.resize(num_points_on_primary_fiber[i]);
    if (record_normMSE)
      fibers[counter].normMSE.resize(num_points_on_primary_fiber[i]);
    fibers[counter].norm.resize(num_points_on_primary_fiber[i]);
    if (record_state)
      fibers[counter].state.resize(num_points_on_primary_fiber[i]);
    if (record_cov)
      fibers[counter].covariance.resize(num_points_on_primary_fiber[i]);

    int k = 0;
    // The first point in the first_half, namely the seed point in the first half, is excluded
    for (int j = static_cast<int>(first_half.position.size()) - 1; j > 0; j--)
    {
      fibers[counter].position[k] = first_half.position[j];
      if (record_rtop1)
        fibers[counter].rtop1[k] = first_half.rtop1[j];
      if (record_rtop2)
        fibers[counter].rtop2[k] = first_half.rtop2[j];
      if (record_rtop3)
        fibers[counter].rtop3[k] = first_half.rtop3[j];
      if (record_rtop_model)
        fibers[counter].rtop_model[k] = first_half.rtop_model[j];
      if (record_rtop_signal)
        fibers[counter].rtop_signal[k] = first_half.rtop_signal[j];
      if (record_free_water)
        fibers[counter].free_water[k] = first_half.free_water[j];
      if (record_w1)
        fibers[counter].w1[k] = first_half.w1[j];
      if (record_w2)
        fibers[counter].w2[k] = first_half.w2[j];
      if (record_w3)
        fibers[counter].w3[k] = first_half.w3[j];
      if (record_w1w2angle)
        fibers[counter].w1w2angle[k] = first_half.w1w2angle[j];
      if (record_w1w3angle)
        fibers[counter].w1w3angle[k] = first_half.w1w3angle[j];
      if (record_Fm1)
        fibers[counter].Fm1[k] = first_half.Fm1[j];
      if (record_lmd1)
        fibers[counter].lmd1[k] = first_half.lmd1[j];
      if (record_Fm2)
        fibers[counter].Fm2[k] = first_half.Fm2[j];
      if (record_lmd2)
        fibers[counter].lmd2[k] = first_half.lmd2[j];
      if (record_Fm3)
        fibers[counter].Fm3[k] = first_half.Fm3[j];
      if (record_lmd3)
        fibers[counter].lmd3[k] = first_half.lmd3[j];
      if (record_varW1)
        fibers[counter].varW1[k] = first_half.varW1[j];
      if (record_varW2)
        fibers[counter].varW2[k] = first_half.varW2[j];
      if (record_varW3)
        fibers[counter].varW3[k] = first_half.varW3[j];
      if (record_varWiso)
        fibers[counter].varWiso[k] = first_half.varWiso[j];
      if (record_normMSE)
        fibers[counter].normMSE[k] = first_half.normMSE[j];
      fibers[counter].norm[k] = first_half.norm[j];
      if (record_state)
        fibers[counter].state[k] = first_half.state[j];
      if (record_cov)
        fibers[counter].covariance[k] = first_half.covariance[j];
      k++;
    }

    for (int j = 0; j < static_cast<int>(second_half.position.size()); j++)
    {
      fibers[counter].position[k] = second_half.position[j];
      if (record_rtop1)
        fibers[counter].rtop1[k] = second_half.rtop1[j];
      if (record_rtop2)
        fibers[counter].rtop2[k] = second_half.rtop2[j];
      if (record_rtop3)
        fibers[counter].rtop3[k] = second_half.rtop3[j];
      if (record_rtop_model)
        fibers[counter].rtop_model[k] = second_half.rtop_model[j];
      if (record_rtop_signal)
        fibers[counter].rtop_signal[k] = second_half.rtop_signal[j];
      if (record_free_water)
        fibers[counter].free_water[k] = second_half.free_water[j];
      if (record_w1)
        fibers[counter].w1[k] = second_half.w1[j];
      if (record_w2)
        fibers[counter].w2[k] = second_half.w2[j];
      if (record_w3)
        fibers[counter].w3[k] = second_half.w3[j];
      if (record_w1w2angle)
        fibers[counter].w1w2angle[k] = second_half.w1w2angle[j];
      if (record_w1w3angle)
        fibers[counter].w1w3angle[k] = second_half.w1w3angle[j];
      if (record_Fm1)
        fibers[counter].Fm1[k] = second_half.Fm1[j];
      if (record_lmd1)
        fibers[counter].lmd1[k] = second_half.lmd1[j];
      if (record_Fm2)
        fibers[counter].Fm2[k] = second_half.Fm2[j];
      if (record_lmd2)
        fibers[counter].lmd2[k] = second_half.lmd2[j];
      if (record_Fm3)
        fibers[counter].Fm3[k] = second_half.Fm3[j];
      if (record_lmd3)
        fibers[counter].lmd3[k] = second_half.lmd3[j];
      if (record_varW1)
        fibers[counter].varW1[k] = second_half.varW1[j];
      if (record_varW2)
        fibers[counter].varW2[k] = second_half.varW2[j];
      if (record_varW3)
        fibers[counter].varW3[k] = second_half.varW3[j];
      if (record_varWiso)
        fibers[counter].varWiso[k] = second_half.varWiso[j];
      if (record_normMSE)
        fibers[counter].normMSE[k] = second_half.normMSE[j];
      fibers[counter].norm[k] = second_half.norm[j];
      if (record_state)
        fibers[counter].state[k] = second_half.state[j];
      if (record_cov)
        fibers[counter].covariance[k] = second_half.covariance[j];
      k++;
    }

    assert(k == num_points_on_primary_fiber[i]);

    counter++;
  }

  assert(counter == num_valid_primary_fibers);
}
