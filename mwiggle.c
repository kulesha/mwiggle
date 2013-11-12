/*
 
Author: Eugene Kulesha 
Contact: kulesh@gmail.com

*/

#include <stdio.h>
#include <stdlib.h>

#include "mw.h"

void usage( void ) {
  printf("Usage:\n\n mwiggle COMMAND FILE PARAMS, e.g\n\nmwiggle create myfile.mw *.wig -t 9060 -a CRGh37 -d \"My study files\"\nmwiggle stats myfile.mw \nmwiggle dump myfile.mw [-r I:1-100] [-t 0,1]\nmwiggle fetch myfile.mw -r I:1-100 [-t 0,1] [-w 20]\n\n");
}

void print_stats(char *fname) {
  META *meta = mw_stats(fname);
  TRACK *tracks = mw_tracks(fname);
  REGION *regions = mw_regions(fname);

  int i = 0; 
  if (meta) {
    printf("*** Stats *** \n");
    while (i < META_NUM) { // maximum number of meta keys
      printf("%s: %s\n", (meta+i)->k, (meta+i)->v);
      i++;
    }    
    free(meta);
  }

  if (regions) {
    REGION *r = regions;
    printf("Regions: ");
    i = 0;
    while (r->size) {
      r = regions + i;
      printf(" %s", r->name);
      i++;
    }
    printf("\n");
    free(regions);
  }

  if (tracks) {
    TRACK *t;
    printf("*** Tracks *** \n");
    i = 0;
    do {
      t = tracks+i;
      printf("Track %d: %s\n", t->id, t->name);
      printf("\t Strand : %c\n", (t->ori > 0) ? '+' : '-');
      printf("\t Values : %f .. %f\n", t->min, t->max);
      printf("\t Desc   : %s\n", t->desc);
      i++;
    } while (t->id);
    free(tracks);
  }


}

int get_data(char *fname, char **opts) {
  char *tracks = NULL;
  char *region = NULL;
  char *window = NULL;
  char *func = NULL;

  int tcount = 0;
  int winsize = 0;

  while (*opts) {
    if ((*opts)[0] == '-') {
      switch((*opts)[1]) {
      case 't': // tracks, e.g -t 0,2
	opts++;
	tracks = *opts;	
	break;
      case 'r': // region, e.g -r I:3-300
	opts++;
	region = *opts;
	break;
      case 'w': // window size in pixels, e.g -w 1000  
	opts++;
	window = *opts;
	winsize = atoi(window);
	if (!winsize) {
	  printf("Error: -w option requires an integer argument, e.g. -w 100\n");
	  return -1;
	}
	break;
     case 'f': // binning function, e.g -f pick|sum 
	opts++;
	func = *opts;
	break;
      default:
	printf("Error: unknown option -%c\n", (*opts)[1]);
	return -1;
      }
    } 
    opts++;
  }

  if (!region) {
    printf("Error: fetch function requires the region parameter, e.g -r I:1-1000\n");
    return -1;
  }

  RESULT *res = mw_fetch(fname, region, tracks, winsize, &tcount);

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
  }
  return 0;
}

int dump_data(char *fname, char **opts) {
  char *tracks = NULL;
  char *region = NULL;

  while (*opts) {
    if ((*opts)[0] == '-') {
      switch((*opts)[1]) {
      case 't': // tracks, e.g -t 0,2
	opts++;
	tracks = *opts;	
	break;
      case 'r': // region, e.g -r I:3-300
	opts++;
	region = *opts;
	break;
      default:
	printf("Error: unknown option -%c\n", (*opts)[1]);
	return -1;
      }
    } 
    opts++;
  }
  mw_dump(fname, region, tracks);
  return 0;
}

int main (int argc, char **argv) {
  char *cmd;
  char *fname;
  
  if (argc < 3) {
    usage();
    return 0;
  }
  argv++;

  cmd = *(argv++);
  fname = *(argv++);

  switch (*cmd) {
  case 'c': // create
    if (mw_create(fname, argv, argc)) {
      usage();
    } else {
      printf("Created %s\n", fname);
      print_stats(fname);
    }
    break;
  case 's': // stats
    print_stats(fname);
    break;
  case 'f': // fetch
    get_data(fname, argv);
    break;
  case 'd': // dump
    dump_data(fname, argv);
    break;
  default:
    printf("Error: Unrecognized command: %s\n", cmd);
  }

  return 0;
}
