#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "sequence.h"
#include "common.h"

// Variant which uses char bases "A","C","G","T" instead of 0-3
// Compute Shannon entropy of seq
//   The natural log is used so the units of entropy are
//   known as "hartley".
//
//   H(X) = -SUM[i-n]( P(x[i])*Logb( P(x[i]) ) )
//
//   The higher the entropy the more disordered the seq.  I.e
//   ACGTACGTACGT is the most disordered 12mer with an
//                entropy of 1.386 hartley
//   AAAAAAAAAGGC is right on the cusp of our cutoff with
//                an entropy of 0.721 hartley
//   AAAAAAAAAAAA and the other monomer 12mer equivs is the most ordered
//                with an entropy of 0 hartley
//
double
compute_seq_entropy(char *seq)
{
  int x;
  int count[4];
  double answer, y;

  int wordLen = strlen(seq);

  for (x = 0; x < 4; x++)
    count[x] = 0;
  for (x = 0; x < wordLen; x++)
    count[(int) char_to_num(seq[x])] += 1;
  answer = 0.0;
  for (x = 0; x < 4; x++)
  {
    if (count[x] == 0)
      continue;
    y = ((double) count[x]) / ((double) wordLen);
    answer += y * log(y);
  }
  // printf("Entropy = %f\n", -answer);
  return -answer;
}


int
rand_int(int n)
{
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do
  {
    rnd = rand();
  }
  while (rnd >= limit);
  return rnd % n;
}
