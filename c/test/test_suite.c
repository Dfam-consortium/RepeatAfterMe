#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "bnw_extend.h"
#include "score_system.h"


int main(int argc, char *argv[]) {
  printf("#\n");
  printf("# Testing bnw_extend...\n");
  printf("#\n");
  bnw_extend_test();
  printf("#\n");
  printf("# Testing score_system..\n");
  printf("#\n");
  score_system_test();
  exit(0);
}
