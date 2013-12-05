/*
 
Author: Eugene Kulesha 
Contact: kulesh@gmail.com

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
      printf("Track ID: %d\n", t->id);
      printf("\t Name : %s\n", t->name);
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

int update_meta(char *fname, char **opts) {
  char *tracks = NULL;
  FILE *f;
  HEADER h;
  int taxid;

  TRACK *ltracks = NULL;
  REGION *regions = NULL;

  int i = 0;

  if ( (f = fopen(fname, "r"))) {
    if (fread(&h, sizeof(HEADER), 1, f)) {
      ltracks =  (TRACK*) malloc(sizeof(TRACK) * h.tracks );
      for (i = 0; i<h.tracks; i++) {
	if ( !fread((ltracks+i), sizeof(TRACK), 1, f)) {
	  printf("Error: could not read from %s\n", fname);
	  fclose(f);
	  free(ltracks);
	  return -1;
	} 
      }

      regions =  (REGION*) malloc(sizeof(REGION) * (h.regions+1) );
      for (i = 0; i<h.regions; i++) {
	if ( !fread((regions+i), sizeof(REGION), 1, f)) {
	  printf("Error: could not read from %s\n", fname);
	  fclose(f);
	  return -1;
	}
      }
      (regions+h.regions)->size = 0; // to indicate the end of the list

    } else {
      printf("Error: could not read %s\n", fname);
      fclose(f);
      return -1;
    }
    
    fclose(f);
  }

  char *name = NULL, *info = NULL;
  UCHAR ori = -1;

  while (*opts) {
    if ((*opts)[0] == '-') {
      switch((*opts)[1]) {
      case 't': // tracks, e.g -t 0,2
	opts++;
	tracks = *opts;	
	break;
      case 'n': // name, e.g -n Study1
	opts++;
	name = *opts;
	break;
      case 's': // strand, e.g -s +
	opts++;
	if (*opts[0] == '-') {
	  ori = 0;
	}  else {
	  ori = 1;
	}
	break;
      case 'i': // strand, e.g -i "this track is for tissue 1"
	opts++;
	info = *opts;
	break;

      case 'd': // description, e.g -d "My first study"
	opts++;
	sprintf(h.desc, *opts);
	break;
      case 'a': // assembly, e.g -a NCBI37
	opts++;
	sprintf(h.assembly, *opts);
	break;
      case 'x': // tax id, e.g -x 9606
	opts++;
	if ((taxid = atoi(*opts)) > 0) {
	  h.taxon = taxid;
	} else {
	  printf("Error: invalid taxonomy id [%s]\n", *opts);
	}
	break;
      default:
	printf("Error: unknown option -%c\n", (*opts)[1]);
	return -1;
      }
    } 
    opts++;
  }


  if (tracks) {
    UCHAR tarray[MAX_TRACKS];
    const char s[2] = ",";
    int tnum, active = -1;
    char *token;

    memset(tarray, 1, MAX_TRACKS);  
    // tracks are specified
    if ((token = strtok(tracks, s))) {
      if ((strchr(token, '-') == NULL)) { // - means all tracks
	memset(tarray, 0, MAX_TRACKS);
	tnum = atoi(token);
	tarray[tnum] = 1;
	active = 1;
	while (token != NULL) {
	  token = strtok(NULL, s);
	  if (token) {
	    tnum = atoi(token);
	    tarray[tnum] = 1;
	    active ++;
	  }
	}
      }
    }

    for (i = 0; i<h.tracks; i++) {
      TRACK *tmp = (ltracks + i);
      if (tarray[ tmp->id ]) {
	if (name) {
	  sprintf(tmp->name, name);
	}
	if (info) {
	  sprintf(tmp->desc, info);
	}
	if (ori >= 0) {
	  tmp->ori = ori;
	}

      }
    }
  }


  if ((f = fopen(fname,"r+"))) {
    TRACK *t;

    fseek(f, 0 , SEEK_SET);
    fwrite(&h, sizeof(HEADER), 1, f);

    if (tracks) {
      for (i = 0; i<h.tracks; i++) {
	t = (ltracks + i);
	fwrite(t, sizeof(TRACK), 1, f);
      } 
    }
    free(ltracks);
    free(regions);
    fclose(f);
  } else {
    printf("Error: failed to open %s for writing\n", fname);
  }



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
  case 'u': // update
    update_meta(fname, argv);
    break;
  default:
    printf("Error: Unrecognized command: %s\n", cmd);
  }

  return 0;
}
