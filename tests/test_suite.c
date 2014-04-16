#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#include "mw.h"


static int test_create() {
  char *fname = "/tmp/test.mw";

  char *args[14] = {
    "data/file1.bedgraph", 
    "data/file2.bedgraph", 
    "data/file3.bedgraph", 
    "data/file4.bedgraph", 
    "data/file5.bedgraph", 
    "data/file6.bedgraph", 
    "-t", "7801", "-a", "AgamP3", "-d", "test", "-v", "5"};

  return mw_create(fname, args, 14);
}

static int test_fetch() {
  char *fname = "/tmp/test.mw";
  int tcount = 0;
  int winsize = 10;

  RESULT *res = mw_fetch(fname, "2L:3040-3050", NULL, winsize, &tcount);

  if (res) {
    int i ;
    RESULT *tr = res;
    while (tcount --) {
      printf("%d: %s\n", tr->t.id, tr->t.name);
      for( i = 0; i < winsize; i++) {
	printf("%f ", tr->v[i]);
      }
      printf("\n");
      free(tr->v);
      tr++;
    }
    free(res);
    return 0;
  }

  return 1;
}



int main() {
  //  assert(false && "My first unit test");

  assert((test_create() == 0) && printf("Test create : Pass\n"));

  assert((test_fetch() == 0) && printf("Test fetch : Pass\n"));
  return 0;
}
