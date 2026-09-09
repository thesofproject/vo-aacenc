/*
 ** Copyright 2003-2010, VisualOn, Inc.
 **
 ** Licensed under the Apache License, Version 2.0 (the "License");
 ** you may not use this file except in compliance with the License.
 ** You may obtain a copy of the License at
 **
 **     http://www.apache.org/licenses/LICENSE-2.0
 **
 ** Unless required by applicable law or agreed to in writing, software
 ** distributed under the License is distributed on an "AS IS" BASIS,
 ** WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ** See the License for the specific language governing permissions and
 ** limitations under the License.
 */
/*******************************************************************************
	File:		band_nrg.c

	Content:	Band/Line energy calculations functions

*******************************************************************************/

#include "basic_op.h"
#include "band_nrg.h"

#ifndef ARMV5E
/********************************************************************************
*
* function name: CalcBandEnergy
* description:   Calc sfb-bandwise mdct-energies for left and right channel
*
**********************************************************************************/
void CalcBandEnergy(const Word32 *mdctSpectrum,
                    const Word16 *bandOffset,
                    const Word16  numBands,
                    Word32       *bandEnergy,
                    Word32       *bandEnergySum)
{
  Word32 i, j;
  Word32 accuSum = 0;

  for (i=0; i<numBands; i++) {
    Word32 accu = 0;
    Word32 start = bandOffset[i];
    Word32 end = bandOffset[i+1];
    j = start;
#if defined(__xtensa__)
    for (; j + 4 <= end; j += 4) {
      Word32 s0 = mdctSpectrum[j+0];
      Word32 s1 = mdctSpectrum[j+1];
      Word32 s2 = mdctSpectrum[j+2];
      Word32 s3 = mdctSpectrum[j+3];
      accu += MULHIGH(s0, s0);
      accu += MULHIGH(s1, s1);
      accu += MULHIGH(s2, s2);
      accu += MULHIGH(s3, s3);
    }
#endif
    for (; j < end; j++)
      accu = L_add(accu, MULHIGH(mdctSpectrum[j], mdctSpectrum[j]));

	accu = L_add(accu, accu);
    accuSum = L_add(accuSum, accu);
    bandEnergy[i] = accu;
  }
  *bandEnergySum = accuSum;
}

/********************************************************************************
*
* function name: CalcBandEnergyMS
* description:   Calc sfb-bandwise mdct-energies for left add or minus right channel
*
**********************************************************************************/
void CalcBandEnergyMS(const Word32 *mdctSpectrumLeft,
                      const Word32 *mdctSpectrumRight,
                      const Word16 *bandOffset,
                      const Word16  numBands,
                      Word32       *bandEnergyMid,
                      Word32       *bandEnergyMidSum,
                      Word32       *bandEnergySide,
                      Word32       *bandEnergySideSum)
{

  Word32 i, j;
  Word32 accuMidSum = 0;
  Word32 accuSideSum = 0;

  for(i=0; i<numBands; i++) {
    Word32 accuMid = 0;
    Word32 accuSide = 0;
    Word32 start = bandOffset[i];
    Word32 end = bandOffset[i+1];
    j = start;
#if defined(__xtensa__)
    for (; j + 4 <= end; j += 4) {
      Word32 l0 = mdctSpectrumLeft[j+0] >> 1, r0 = mdctSpectrumRight[j+0] >> 1;
      Word32 l1 = mdctSpectrumLeft[j+1] >> 1, r1 = mdctSpectrumRight[j+1] >> 1;
      Word32 l2 = mdctSpectrumLeft[j+2] >> 1, r2 = mdctSpectrumRight[j+2] >> 1;
      Word32 l3 = mdctSpectrumLeft[j+3] >> 1, r3 = mdctSpectrumRight[j+3] >> 1;
      Word32 m0 = l0 + r0, s0 = l0 - r0;
      Word32 m1 = l1 + r1, s1 = l1 - r1;
      Word32 m2 = l2 + r2, s2 = l2 - r2;
      Word32 m3 = l3 + r3, s3 = l3 - r3;
      accuMid += MULHIGH(m0, m0);
      accuSide += MULHIGH(s0, s0);
      accuMid += MULHIGH(m1, m1);
      accuSide += MULHIGH(s1, s1);
      accuMid += MULHIGH(m2, m2);
      accuSide += MULHIGH(s2, s2);
      accuMid += MULHIGH(m3, m3);
      accuSide += MULHIGH(s3, s3);
    }
#endif
    for (; j < end; j++) {
      Word32 specm, specs;
      Word32 l, r;

      l = mdctSpectrumLeft[j] >> 1;
      r = mdctSpectrumRight[j] >> 1;
      specm = l + r;
      specs = l - r;
      accuMid = L_add(accuMid, MULHIGH(specm, specm));
      accuSide = L_add(accuSide, MULHIGH(specs, specs));
    }

	accuMid = L_add(accuMid, accuMid);
	accuSide = L_add(accuSide, accuSide);
	bandEnergyMid[i] = accuMid;
    accuMidSum = L_add(accuMidSum, accuMid);
    bandEnergySide[i] = accuSide;
    accuSideSum = L_add(accuSideSum, accuSide);

  }
  *bandEnergyMidSum = accuMidSum;
  *bandEnergySideSum = accuSideSum;
}

#endif
