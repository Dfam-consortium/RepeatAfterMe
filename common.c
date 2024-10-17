#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "sequence.h"
#include "common.h"


void
print_small_kmer_int(uint32_t val, int kmer_k)
{
  int i;
  int mask = 3;
//printf("Why is this: %d\n", (int)pow((double)4,(double)7)*3);
  if ( kmer_k > 1 )
    mask = mask * (int)(pow((double)4,(double)(kmer_k-1)));
  for (i = kmer_k-1; i >= 0; i-- )
  {
    //printf("%c %d mask = %d", num_to_char((char)((val & mask)>>(i*2))), val, mask);
    printf("%c", num_to_char((char)((val & mask)>>(i*2))));
    mask = mask >> 2;
  }
}

// TODO: move this to hash_functions.c
// small_kmer_to_int
//
//   Convert a kmer for k's between 1 and 15 into
//   a 32bit integer.  Bases are packed as 2 bits
//   with the first base in the MSB using the
//   function ['A'=0,'C'=1,'G'=2','T'=3] to map
//   bases to bits.
//   If the user specifies a k > 15 or there are
//   N's in the sequence this function returns
//   the maximum 32bit integer (0xFFFFFFFF) as
//   an error code.  If rc is set reverse complement
//   the sequence.
//
uint32_t
small_kmer_to_int (char *kmer, int kmer_k, int rc)
{
  int i;
  uint32_t kmer_int = 0;

  if ( kmer_k < 16 )
  {
    if ( rc )
    {
      for (i = kmer_k-1; i >= 0; i-- )
      {
        kmer_int = kmer_int << 2;
        // A=00 C=01,G=10,T=11
        // a=0100 c=0101, g=0110, t=0111
        //  xor with 0b11
        // T=11 G=10,C=01,A=00
        if ( kmer[i] <= 7 )
          kmer_int += (kmer[i]&3)^3; // xor 0b11
        else
          return((uint32_t)0xFFFFFFFF);
      }
    }else{
      for (i = 0; i < kmer_k; i++ )
      {
        kmer_int = kmer_int << 2;
        if ( kmer[i] <= 7 )
          kmer_int += (kmer[i]&3);
        else
          return((uint32_t)0xFFFFFFFF);
      }
    }
    return(kmer_int);
  }else {
    return((uint32_t)0xFFFFFFFF);
  }
}

// small_kmer_seq_to_int
//
//   Convert a kmer for k's between 1 and 15 into
//   a 32bit integer.  Bases are packed as 2 bits
//   with the first base in the MSB using the
//   function ['A'=0,'C'=1,'G'=2','T'=3] to map
//   bases to bits.
//   If the user specifies a k > 15 or there are
//   N's in the sequence this function returns
//   the maximum 32bit integer (0xFFFFFFFF) as
//   an error code.  If rc is set reverse complement
//   the sequence.
//
uint32_t
small_kmer_seq_to_int (char *kmer, int kmer_k, int rc)
{
  int i;
  uint32_t kmer_int = 0;
  char ch;

  if ( kmer_k < 16 )
  {
    if ( rc )
    {
      for (i = kmer_k-1; i >= 0; i-- )
      {
	kmer_int = kmer_int << 2;
	// A=00 C=01,G=10,T=11
	// a=0100 c=0101, g=0110, t=0111
	//  xor with 0b11
	// T=11 G=10,C=01,A=00
        ch = char_to_num(kmer[i]);
	if ( ch < 8 )
	  kmer_int += (ch&3)^3; // xor 0b11
	else
	  return((uint32_t)0xFFFFFFFF);
      }
    }else{
      for (i = 0; i < kmer_k; i++ )
      {
	kmer_int = kmer_int << 2;
        ch = char_to_num(kmer[i]);
	if ( ch < 8 )
	  kmer_int += (ch&3);
	else
	  return((uint32_t)0xFFFFFFFF);
      }
    }
    return(kmer_int);
  }else {
    return((uint32_t)0xFFFFFFFF);
  }
}

// large_kmer_to_int
//
//   Convert a kmer for k's between 1 and 31 into
//   a 64bit integer.  Bases are packed as 2 bits
//   with the first base in the MSB using the
//   function ['A'=0,'C'=1,'G'=2','T'=3] to map
//   bases to bits.
//   If the user specifies a k > 31 or there are
//   N's in the sequence this function returns
//   the maximum 64bit integer (0xFFFFFFFFFFFFFFFF) as
//   an error code.  If rc is set, reverse complement
//   the sequence.
//
uint64_t
large_kmer_to_int (char *kmer, int kmer_k, int rc)
{
  int i;
  uint32_t kmer_int = 0;

  if ( kmer_k < 32 )
  {
    if ( rc )
    {
      for (i = kmer_k-1; i >= 0; i-- )
      {
        kmer_int = kmer_int << 2;
        // A=00 C=01,G=10,T=11
        // a=0100 c=0101, g=0110, t=0111
        //  xor with 0b11
        // T=11 G=10,C=01,A=00
        if ( kmer[i] <= 7 )
          kmer_int += (kmer[i]&3)^3; // xor 0b11
        else
          return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
    }else{
      for (i = 0; i < kmer_k; i++ )
      {
        kmer_int = kmer_int << 2;
        if ( kmer[i] <= 7 )
          kmer_int += (kmer[i]&3);
        else
          return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
    }
    return(kmer_int);
  }else {
    return((uint64_t)0xFFFFFFFFFFFFFFFF);
  }
}

//
// large_kmer_pair_to_int
//
//   Convert a kmer for k's between 1 and 31 into
//   a 64bit integer.  Bases are packed as 2 bits
//   with the first base in the MSB using the
//   function ['A'=0,'C'=1,'G'=2','T'=3] to map
//   bases to bits.
//   If the user specifies a k > 31 or there are
//   N's in the sequence this function returns
//   the maximum 64bit integer (0xFFFFFFFFFFFFFFFF) as
//   an error code.  If rc is set, reverse complement
//   the sequence.
//
uint64_t
large_kmer_pair_to_int (char *leftKmer, char *rightKmer, int kmerSize, int rc)
{
  int i;
  uint32_t kmer_int = 0;

  if ( kmerSize < 32 )
  {
    if ( rc )
    {
      for (i = kmerSize-1; i >= 0; i-- )
      {
        kmer_int = kmer_int << 2;
        // A=00 C=01,G=10,T=11
        // a=0100 c=0101, g=0110, t=0111
        //  xor with 0b11
        // T=11 G=10,C=01,A=00
        if ( rightKmer[i] <= 7 )
          kmer_int += (rightKmer[i]&3)^3; // xor 0b11
        else
          return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
      for (i = kmerSize-1; i >= 0; i-- )
      {
        kmer_int = kmer_int << 2;
        // A=00 C=01,G=10,T=11
        // a=0100 c=0101, g=0110, t=0111
        //  xor with 0b11
        // T=11 G=10,C=01,A=00
        if ( leftKmer[i] <= 7 )
          kmer_int += (leftKmer[i]&3)^3; // xor 0b11
        else
          return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
    }else{
      for (i = 0; i < kmerSize; i++ )
      {
        kmer_int = kmer_int << 2;
        if ( leftKmer[i] <= 7 )
          kmer_int += (leftKmer[i]&3);
        else
          return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
      for (i = 0; i < kmerSize; i++ )
      {
        kmer_int = kmer_int << 2;
        if ( rightKmer[i] <= 7 )
          kmer_int += (rightKmer[i]&3);
        else
          return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
    }
    return(kmer_int);
  }else {
    return((uint64_t)0xFFFFFFFFFFFFFFFF);
  }
}

// large_kmer_seq_to_int
//
//   Convert a kmer for k's between 1 and 31 into
//   a 64bit integer.  Bases are packed as 2 bits
//   with the first base in the MSB using the
//   function ['A'=0,'C'=1,'G'=2','T'=3] to map
//   bases to bits.
//   If the user specifies a k > 31 or there are
//   N's in the sequence this function returns
//   the maximum 64bit integer (0xFFFFFFFFFFFFFFFF) as
//   an error code.  If rc is set, reverse complement
//   the sequence.
//
uint64_t
large_kmer_seq_to_int (char *kmer, int kmer_k, int rc)
{
  int i;
  uint32_t kmer_int = 0;
  char ch;

  if ( kmer_k < 32 )
  {
    if ( rc )
    {
      for (i = kmer_k-1; i >= 0; i-- )
      {
	kmer_int = kmer_int << 2;
	// A=00 C=01,G=10,T=11
	// a=0100 c=0101, g=0110, t=0111
	//  xor with 0b11
	// T=11 G=10,C=01,A=00
        ch = char_to_num(kmer[i]);
	if ( ch < 8 )
	  kmer_int += (ch&3)^3; // xor 0b11
	else
	  return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
    }else{
      for (i = 0; i < kmer_k; i++ )
      {
	kmer_int = kmer_int << 2;
        ch = char_to_num(kmer[i]);
	if ( ch < 8 )
	  kmer_int += (ch&3);
	else
	  return((uint64_t)0xFFFFFFFFFFFFFFFF);
      }
    }
    return(kmer_int);
  }else {
    return((uint64_t)0xFFFFFFFFFFFFFFFF);
  }
}



// Compute Shannon entropy of kmer
//   The natural log is used so the units of entropy are
//   known as "hartley".
//
//   H(X) = -SUM[i-n]( P(x[i])*Logb( P(x[i]) ) )
//
//   The higher the entropy the more disordered the kmer.  I.e
//   ACGTACGTACGT is the most disordered 12mer with an
//                entropy of 1.386 hartley
//   AAAAAAAAAGGC is right on the cusp of our cutoff with
//                an entropy of 0.721 hartley
//   AAAAAAAAAAAA and the other monomer 12mer equivs is the most ordered
//                with an entropy of 0 hartley
//
//  NOTE: returns -1.0 if there is an error with conversion
//
double
compute_entropy(char *kmer, int kmer_k)
{
  int x;
  int count[4];
  double answer, y;

  for (x = 0; x < 4; x++)
    count[x] = 0;
  // ACGT=0-3, acgt=4-7
  // and 3 = a->A, c->C, g->G, t->T
  // N's cause an error
  for (x = 0; x < kmer_k; x++)
    if ( kmer[x] < 8 )
      count[kmer[x]&3] += 1;
    else
      return( -1.0 );
  answer = 0.0;
  for (x = 0; x < 4; x++)
  {
    if (count[x] == 0)
      continue;
    y = ((double) count[x]) / ((double) kmer_k);
    answer += y * log(y);
  }
  return -answer;
}

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
//  NOTE: returns -1.0 if there is an error with conversion
double
compute_seq_entropy(char *seq)
{
  int x;
  int count[4];
  double answer, y;

  int wordLen = strlen(seq);
  char ch;

  for (x = 0; x < 4; x++)
    count[x] = 0;

  for (x = 0; x < wordLen; x++) {
    ch = char_to_num(seq[x]);
    if ( ch < 8 )
      count[ch&3] += 1;
    else
      return(-1.0);
  }

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
