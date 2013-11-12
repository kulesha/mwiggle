/*
 
Author: Eugene Kulesha 
Contact: kulesh@gmail.com

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>

#include "mw.h"

#define MAX_PATH 1024

TRECORD *troot = NULL, *tcur = NULL;
RRECORD *rroot = NULL, *rcur = NULL;

int create_index( char *fname, UINT taxon, char *assembly, char *desc ) {
  HEADER h;
  char hfile[250];
  FILE *hf;
  char cmd[1024];
  double minV, maxV;

  sprintf(hfile, "%s.head", "test01");

  sprintf(h.format, FORMAT_MWIG);
  h.v_major = VERSION_MAJOR;
  h.v_minor = VERSION_MINOR;
  h.taxon = taxon;
  sprintf(h.assembly, assembly);

  sprintf(h.desc, desc);

  h.tracks = 0;
  tcur = troot;
  if (tcur) {
    minV = tcur->data.min;
    maxV = tcur->data.max;
  }
  while (tcur) {
    if (minV > tcur->data.min) {
      minV = tcur->data.min;
    }

    if (maxV < tcur->data.max) {
      maxV = tcur->data.max;
    }
      
    tcur = tcur->next;
    h.tracks++;
  }

  h.regions = 0;
  h.min = minV;
  h.max = maxV;

  rcur = rroot;
  while (rcur) {
    rcur = rcur->next;
    h.regions++;
  }

  h.offset = sizeof(HEADER) + sizeof(TRACK) * h.tracks + sizeof(REGION) * h.regions;
  //  printf("\nFormat: %s\nVersion: %d.%d\nTracks: * %d\nRegions: %d\nOffset: %ld\nTaxon ID: %d\nAssembly: %s\nValues: %f .. %f\n\n", h.format, h.v_major, h.v_minor, h.tracks, h.regions, h.offset, h.taxon, h.assembly, h.min, h.max);


  if ((hf = fopen(hfile,"w+"))) {
    fwrite(&h, sizeof(HEADER), 1, hf);

    tcur = troot;
    while (tcur) {
      fwrite(tcur, sizeof(TRACK), 1, hf);
      //      printf("%d: %s : %s ( %f .. %f )\n", tcur->data.id, tcur->data.name, tcur->data.desc, tcur->data.min, tcur->data.max);
      tcur = tcur->next;
      free(troot);
      troot = tcur;
    }


    rcur = rroot;
    while (rcur) {
      fwrite(rcur, sizeof(REGION), 1, hf);
      //      printf("%s : %ld : %ld\n", rcur->data.name, rcur->data.offset, rcur->data.size);
      rcur = rcur->next;
      free(rroot);
      rroot = rcur;
    }
    fclose(hf);
  }

  sprintf(cmd, "cat test01.head test01.body > %s", fname);
  //  printf(cmd);
  //printf("\n");
  system(cmd);


  sprintf(cmd, "rm *.chr* test01.* ");
  //  printf(cmd);
  //printf("\n");
  system(cmd);
  
  return 0;
}

int merge_regions ( int total ) {
  struct dirent *e;
  DIR *d = opendir(".");
  FILE *r, *b;
  char bfile[250];
  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  REGION a;

  double* Values = malloc(sizeof(double) * total);
  printf("Merging regions ...\n");
  if (d) {
    sprintf(bfile, "%s.body", "test01");
    if ((b = fopen(bfile,"w+"))) {
      while ( (e = readdir(d)) != NULL) {
	if (strstr(e->d_name, "_sorted")) {
	  //	    sprintf(a.name, e->d_name);
	  char *pchr = strrchr(e->d_name, '.');
	  int l = pchr - e->d_name;
	  strncpy(a.name, e->d_name, l);
	  a.name[l] = '\0';

	  printf(" %s ( %s)\n", e->d_name, a.name);
	  if ( (r = fopen(e->d_name, "r"))) {
	    ULONG pos, curpos = -1;
	    float val;
	    int ind;
	    int num;

	    rcur = (RRECORD *) malloc(sizeof(RRECORD));
	    a.offset = ftell(b);
	    
	    while((read = getline(&line, &len, r)) > 0 ) {
	      if ((num = sscanf(line, "%ld %f %d", &pos, &val, &ind) == 3)) {
		if ( curpos != pos ) {
		  if (curpos != -1) {
		    fwrite(&curpos, sizeof(curpos), 1, b);
		    fwrite(Values, sizeof(double), total, b);
		  }
		  curpos = pos;
		  memset(Values, 0, total);
		} 
		*(Values + ind) = val;
	      }
	    }

	    a.size = ftell(b)-a.offset;
	    memcpy(&rcur->data, &a, sizeof(REGION));

	    rcur->next = rroot;
	    rcur->prev = NULL;
	    
	    if (rroot) {
	      rroot->prev = rcur;
	    }
	    rroot = rcur;

	    fclose(r);
	  }
	}
      }
      fclose(b);
    }
    
  }

  if (Values) {
    free(Values);
  }

  return 0;
}

int sort_regions ( int total ) {
  
  struct dirent *e;
  DIR *d = opendir(".");
  char cmd[200];

  printf("Sorting regions ... \n");
  if (d) {
    while ( (e = readdir(d)) != NULL) {
      if (strstr(e->d_name, ".chr")) {
	printf(" %s\n", e->d_name);
	sprintf(cmd, "sort -g %s > %s_sorted", e->d_name, e->d_name);
	//	printf("%s\n", cmd);
	system(cmd);
      }
    }
    closedir(d);
  }
  return 0;
}


int read_wig (char *fname, int idx, int total) {
  FILE *src;
  FILE *header = NULL;

  char hfile[255];
  char bfile[255];

  ssize_t read;
  size_t len = 0;
  char *line = NULL;
  
  sprintf(hfile, "%s.head", fname);
  sprintf(bfile, "%s.body", fname);
  int state = 0; // reading the step line
  char strStep[50], strChr[50], strSpan[50];
  int step = 0;
  int num ;
  long pos;
  float val;
  char region[50];
  char curRegion[50] = "undefined";
  char description[255] = "";
  TRACK a;
  double minV;
  double maxV;
  int vSet = 0;
  char *slash, *dot;

  if ((src = fopen(fname, "r"))) {
    printf(" * adding %d of %d (%s) \n", idx+1, total, fname);

    while((read = getline(&line, &len, src)) > 0 ) {
      int f = 1;
      while (f) {
	switch(state) {
	case 0:
	  if ((num = sscanf(line, "%s %s %s", strStep, strChr, strSpan) == 3)) {
	    sscanf(strChr, "chrom=%s", region);
	    //	    printf("A: %s  * %s * %s \n ", strChr, region, curRegion);
	    if (strstr(strStep, "fixed") != NULL) {
	      step = 1;
	      state = 1; // fixed step
	    } else {
	      state = 2; // variable step
	    }
	    f = 0;
	    //	    printf("B: %s * %s \n ", region, curRegion);
	    if (strcmp(curRegion, region) != 0) {
	      //	      	      printf(" write %s \n", region);
	      strcpy(curRegion, region);
	      sprintf(hfile, "%s.chr", region);
	      if (header) {
		fclose(header);
	      }

	      if ((header = fopen(hfile, "a+"))) {
		//		printf("opened %s (%d) \n", hfile, state);
	      } else {
		printf(" Cannot open %s ", hfile);
		exit(-1);
	      }
	    }
	  }
	  break;
	case 1:
	  f = 0;
	  break;
	case 2:
	  if ((num = sscanf(line, "%ld %f", &pos, &val) == 2)) {
	    f = 0;
	    if (vSet) {
	      if (minV > val ) {
		minV = val ;
	      }
	      if (maxV < val ) {
		maxV = val ;
	      }
	    } else {
	      vSet = 1;
	      minV = val;
	      maxV = val;
	    }
	    fprintf(header, "%ld %f %d\n", pos, val, idx);
	  } else {
	    state = 0;
	  }
	  break;
	default:
	  f = 0;
	  
	}
      }
    }
    fclose(header);
    fclose(src);

    tcur = (TRECORD *) malloc(sizeof(TRECORD));
    slash = strrchr(fname, '/');
    
    if (slash) {
      sprintf(a.desc, slash+1);
    } else {
      sprintf(a.desc, fname);
    }
    dot = strrchr(a.desc, '.');
    if (dot) {
      strncpy(a.name, a.desc, (dot - a.desc));
    } else {
      sprintf(a.name, a.desc);
    }
    sprintf(a.desc, description);
    a.id = idx;
    a.min = minV;
    a.max = maxV;
    a.ori = 1;
    memcpy(&tcur->data, &a, sizeof(TRACK));
    tcur->next = troot;
    tcur->prev = NULL;

    if (troot) {
      troot->prev = tcur;
    }
    troot = tcur;
  }

  if (line) {
    free(line);
  }
  printf("\t values range: %f .. %f\n", minV, maxV);
  return 0;
}


int mw_create(char *fname, char **argv, int argc) {
  int taxon = 0;
  char *tmp = NULL;
  char *assembly = NULL;
  char *desc = NULL;
  char files[MAX_TRACKS][MAX_PATH];

  int i, num = 0;
  
  while (*argv) {
    if ((*argv)[0] == '-') {
      switch((*argv)[1]) {
      case 't':
	argv++;
	tmp = *argv;
	if ((taxon = atoi(tmp)) == 0) {
	  printf("Error: invalid Taxonomy ID %s\n", tmp);
	}
	
	break;
      case 'a':
	argv++;
	assembly = *argv;
	break;
      case 'd':
	argv++;
	desc = *argv;
	break;
      default:
	printf("Error: unknown option -%c\n", (*argv)[1]);
	return -1;
      }
      argv++;
    } else {
       strcpy(files[num], *argv++);
       num ++;
    }
  }
  
  if (!taxon || !assembly || !desc) {
    printf("Error: Taxonomy ID, Assembly and Description are required to create a MW file\n");
    return -1;
  }

  printf("Taxon : %d\n", taxon);
  printf("Assembly : %s\n", assembly);
  printf("Description: %s\n", desc);

  for(i= 0; i < num ; i++) {
    read_wig(files[i], i, num);
  }

  sort_regions(num);
  merge_regions(num);
  create_index(fname, taxon, assembly, desc);

  return 0;
}

META *mw_stats(char *fname) {
  FILE *f;
  HEADER h;
  META *stats =  malloc(sizeof(META) * META_NUM );
  int i = 0;

  if ( (f = fopen(fname, "r"))) {
    if (fread(&h, sizeof(HEADER), 1, f)) {
      sprintf(stats[i].k, "Version.Major"); sprintf(stats[i++].v, "%d", h.v_major);
      sprintf(stats[i].k, "Version.Minor"); sprintf(stats[i++].v, "%d", h.v_minor);
      sprintf(stats[i].k, "TaxID"); sprintf(stats[i++].v, "%d", h.taxon);
      sprintf(stats[i].k, "Assembly");sprintf(stats[i++].v,"%s", h.assembly);
      sprintf(stats[i].k, "Value.Min");sprintf(stats[i++].v, "%f", h.min);
      sprintf(stats[i].k, "Value.Max");sprintf(stats[i++].v, "%f", h.max);
      sprintf(stats[i].k, "Tracks count");sprintf(stats[i++].v, "%d", h.tracks);
      sprintf(stats[i].k, "Regions count");sprintf(stats[i++].v, "%d", h.regions);
      sprintf(stats[i].k, "Desc");sprintf(stats[i++].v, "%s", h.desc);
      sprintf(stats[i].k, "Format"); strcpy(stats[i++].v,  h.format);
    } else {
      printf("Error: could not read %s\n", fname);
      fclose(f);
      return NULL;
    }
    
    fclose(f);
    return stats;
  }
  printf("Error: could not open %s\n", fname);
  return NULL;
}

TRACK *mw_tracks (char *fname) {
  FILE *f;
  HEADER h;
  TRACK *tracks = NULL;
  int i = 0;

  if ( (f = fopen(fname, "r"))) {
    if (fread(&h, sizeof(HEADER), 1, f)) {
      tracks =  malloc(sizeof(TRACK) * h.tracks );
      for (i = 0; i<h.tracks; i++) {
	if ( !fread((tracks+i), sizeof(TRACK), 1, f)) {
	  printf("Error: could not read from %s\n", fname);
	  fclose(f);
	  free(tracks);
	  return NULL;
	} 
      }
    } else {
      printf("Error: could not read %s\n", fname);
      fclose(f);
      return NULL;
    }
    
    fclose(f);
    return tracks;
  }
  printf("Error: could not open %s\n", fname);
  return tracks;
}

REGION *mw_regions (char *fname) {
  FILE *f;
  HEADER h;
  TRACK t;
  REGION *regions = NULL;
  int i = 0;

  if ( (f = fopen(fname, "r"))) {
    if (fread(&h, sizeof(HEADER), 1, f)) {
      for (i = 0; i<h.tracks; i++) {
	if ( !fread(&t, sizeof(TRACK), 1, f)) {
	  printf("Error: could not read from %s\n", fname);
	  fclose(f);
	  return NULL;
	}
      }
      regions =  malloc(sizeof(REGION) * (h.regions+1) );
      for (i = 0; i<h.regions; i++) {
	if ( !fread((regions+i), sizeof(REGION), 1, f)) {
	  printf("Error: could not read from %s\n", fname);
	  fclose(f);
	  return NULL;
	}
      }
      (regions+h.regions)->size = 0; // to indicate the end of the list

    } else {
      printf("Error: could not read %s\n", fname);
      fclose(f);
      return NULL;
    }
    
    fclose(f);
    return regions;
  }
  printf("Error: could not open %s\n", fname);
  return NULL;
}


RESULT *mw_fetch(char *fname, char *region, char *tracks, int winsize, int *tcount) {
  FILE *f;
  HEADER h;
  REGION r;

  ULONG data_offset;
  ULONG region_offset = -1;
  ULONG region_size = -1;

  char region_name[255] = "", region_start[255];
  ULONG start;
  ULONG end;
  ULONG pos = 0;

  size_t cpos = 0;


  UCHAR tarray[MAX_TRACKS];
  const char s[2] = ",";
  int tnum, active = -1;
  char *token;

  RESULT *res = NULL, *tr;
  int i, j;

  char bfunc = 'p'; // p - pick , s - sum, a - avg ? 

  if (region) {
    char *colon = strchr(region, ':');
    char *hiphen = strchr(region, '-');
    memset(region_name, 0, 255);
    if (colon) {
      strncpy(region_name, region, (colon-region));  
      strncpy(region_start, colon+1, hiphen-colon);
      start = atol(region_start);
      end = atol(hiphen+1);
      if (!start || !end) {
	printf("Error: wrong region format %s\n", region);
	return NULL;
      }
    } else {
      sprintf(region_name, region);
    }
    
    //    printf(" Region %s:%ld-%ld\n", region_name, start, end);
  } else{
    printf("Error: mw_fetch requires a region parameter, e.g. I:1-100\n");
    return NULL;
  }

  // assume all tracks have to be read
  memset(tarray, 1, MAX_TRACKS);  
  if (tracks) {
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
  }

  if ( (f = fopen(fname, "r"))) {
    if (fread(&h, sizeof(HEADER), 1, f)) {
      data_offset = h.offset;
      if (active < 0) { // all tracks 
	active = h.tracks;
      }
      *tcount = active;

      res = (RESULT*) malloc (sizeof (RESULT) * (active + 1));
      j = 0;

      for (i = 0; i<h.tracks; i++) {
	tr = (res + j);
	if (fread(&tr->t, sizeof(TRACK), 1, f)) {
	  if (tarray[i]) {
	    j++;
	  }
	} else {
	  printf("Error: could not read %s\n", fname);
	}
      }

      for (i = 0; i<h.regions; i++) {
	if (fread(&r, sizeof(REGION), 1, f)) {
	  //	  printf(" %s : %ld\n", r.name, r.offset);
	  if (strcmp(r.name, region_name) == 0) {
	    region_offset = r.offset + data_offset;
	    region_size = r.size + region_offset;
	  }
	} else {
	  printf("Error: could not read %s\n", fname);
	}
      }

      if (region_offset != -1) {
	double *values = malloc(sizeof(double) * (h.tracks));
	fseek(f, region_offset, SEEK_SET);
	cpos = ftell(f);
	  
	while (!feof(f) && (cpos < region_size) && (pos < start)) {
	  if (fread(&pos, sizeof (ULONG) , 1, f)) {
	  }
	  if (fread(values, sizeof (double) , h.tracks, f)) {
	  }
	  cpos = ftell(f);
	}

	if (pos > start) {
	  if (winsize) {
	    float step = (end - start + 1) / (float)winsize;
	    int v = 0;
	    int idx;

	    double *binvalues = malloc( sizeof(double) * winsize * h.tracks);
	    memset(binvalues, 0, sizeof(double) * winsize * h.tracks);
	    while (!feof(f) && (cpos < region_size) && (pos < end)) {
	      v ++;
	      idx = (pos - start) / step;
	      for (i = 0; i < h.tracks; i++) {
		if (tarray[i]) {
		  double *cv = binvalues + winsize * i + idx;
		  switch (bfunc) {
		  case 'a': // aggregate
		    values[i] += *cv;
		    break;
		  default: // pick
		    if (fabs(values[i]) > fabs(*cv)) {
		      *cv = values[i];
		    }
		  }
		}
	      }
	      if (fread(&pos, sizeof (ULONG) , 1, f)) {
	      }
	      if (fread(values, sizeof (double) , h.tracks, f)) {
	      }
	      cpos = ftell(f);
	    }

	    i = 0;
	    for (j = 0; j < h.tracks; j++) {
	      if (tarray[j]) {
		RESULT* ptr = (res + i);
		ptr->v = (double*) malloc(sizeof(double) * winsize);
		memcpy(ptr->v, (binvalues + j*winsize), winsize * sizeof(double));
		i ++;
	      }
	    }	      	

	    
	    return res;
	  }

	  while (!feof(f) && (cpos < region_size) && (pos < end)) {
	    printf("%ld", pos);
	    for (i = 0; i < h.tracks; i++) {
	      if (tarray[i]) {
		printf(" %f", values[i]);
	      }
	    }
	    printf("\n");

	    if (fread(&pos, sizeof (ULONG) , 1, f)) {}
	    if (fread(values, sizeof (double) , h.tracks, f)) {}
	    cpos = ftell(f);
	  }
	}

	free(values);
      } else {
	printf("Error: could not find the region %s\n", region_name);
      }
      free(res);
	
    } else {
      printf("Error: could not read %s\n", fname);
    }
    fclose(f);
  }

  return NULL;
}

int mw_dump(char *fname, char *region, char *tracks) {
  FILE *f;
  HEADER h;
  REGION r;
  TRACK trackarray [MAX_TRACKS];

  ULONG data_offset;
  ULONG region_offset = -1;
  ULONG region_end = -1;

  char region_name[255] = "", region_start[255];
  ULONG start = 1;
  ULONG end = -1;
  ULONG pos = 0;

  size_t cpos = 0;


  UCHAR tarray[MAX_TRACKS];
  const char s[2] = ",";
  int tnum, active = -1;
  char *token;


  int i, j;

  FILE *farray[MAX_TRACKS];

  // assume all tracks have to be read
  memset(tarray, 1, MAX_TRACKS);  
  if (tracks) {
    printf(" Dumping tracks %s\n", tracks);
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
  }

  memset(region_name, 0, 255);

  if (region) {
    char *colon = strchr(region, ':');
    char *hiphen = strchr(region, '-');

    if (colon) {
      strncpy(region_name, region, (colon-region));  
      strncpy(region_start, colon+1, hiphen-colon);
      start = atol(region_start);
      end = atol(hiphen+1);
      if (!start || !end) {
	printf("Error: wrong region format %s\n", region);
	return -1;
      }
    } else {
      sprintf(region_name, region);
    }
    
    printf(" Region %s:%ld-%ld\n", region_name, start, end);
  }


  if ( (f = fopen(fname, "r"))) {
    if (fread(&h, sizeof(HEADER), 1, f)) {
      data_offset = h.offset;
      if (active < 0) { // all tracks 
	active = h.tracks;
      }

      j = 0;

      for (i = 0; i<h.tracks; i++) {
	if (fread(&trackarray[i], sizeof(TRACK), 1, f)) {
	  if (tarray[i]) {
	    char dname[1024];
	    sprintf(dname, "%s-%d.%s", trackarray[i].name, trackarray[i].ori, h.format);

	    if ((farray[i] = fopen(dname, "w"))) {
	    } else {
	      printf("Error: could not write into %s\n", dname);
	    }
	  }
	} else {
	  printf("Error: could not read %s\n", fname);
	}
      }

      if (region) {
	for (i = 0; i<h.regions; i++) {
	  if (fread(&r, sizeof(REGION), 1, f)) {
	    if (strcmp(r.name, region_name) == 0) {
	      printf(" Found %s. %ld  / %ld\n", r.name, r.offset, r.size);
	      region_offset = r.offset + data_offset;
	      region_end = r.size + region_offset;
	    }
	  } else {
	    printf("Error: could not read %s\n", fname);
	  }
	}

	if (region_offset != -1) {
	  double *values;
	  
	  fseek(f, region_offset, SEEK_SET);
	  cpos = ftell(f);

	  values = malloc(sizeof(double) * (h.tracks));

	  while (!feof(f) && (cpos < region_end) && (pos < start)) {
	    if (fread(&pos, sizeof (ULONG) , 1, f)) {
	    }
	    if (fread(values, sizeof (double) , h.tracks, f)) {
	    }
	    cpos = ftell(f);
	  }
	  
	  if (pos > start) {
	    for (i = 0; i < h.tracks; i++) {
	      if (tarray[i]) {
		fprintf(farray[i], "variableStep chrom=%s span=1\n", region_name);
	      }
	    }
	    
	    
	    while (!feof(f) && (cpos < region_end) && (pos < end)) {
	      for (i = 0; i < h.tracks; i++) {
		if (tarray[i]) {
		  fprintf(farray[i], "%ld %f\n", pos, values[i]);
		}
	      }
	      if (fread(&pos, sizeof (ULONG) , 1, f)) {}
	      if (fread(values, sizeof (double) , h.tracks, f)) {}
	      cpos = ftell(f);
	    }
	  }

	  free(values);

	} else {
	  printf("Error: could not find the region %s\n", region);
	  return -1;
	}
      } else { // all regions
	double *values;
	ULONG tpos;

	values = malloc(sizeof(double) * (h.tracks));

	region_offset = sizeof(HEADER) + sizeof(TRACK) * h.tracks;
	fseek(f, region_offset, SEEK_SET);
	tpos = ftell(f);
	for (j = 0; j<h.regions; j++) {
	  fseek(f, tpos, SEEK_SET);
	  if (fread(&r, sizeof(REGION), 1, f)) {
	    tpos = ftell(f);
	    //	    printf(" Found %s. %ld  / %ld\n", r.name, r.offset, r.size);
	    region_offset = r.offset + data_offset;
	    region_end = r.size + region_offset;

	    fseek(f, region_offset, SEEK_SET);
	    cpos = ftell(f);

	    while (!feof(f) && (cpos < region_end)) {
	      for (i = 0; i < h.tracks; i++) {
		if (tarray[i]) {
		  fprintf(farray[i], "%ld %f\n", pos, values[i]);
		}
	      }
	      if (fread(&pos, sizeof (ULONG) , 1, f)) {}
	      if (fread(values, sizeof (double) , h.tracks, f)) {}
	      cpos = ftell(f);
	    }


	  } else {
	    printf("Error: could not read %s\n", fname);
	  }
	}
	    free(values);
      
      }


      for (i = 0; i<h.tracks; i++) {
	if (farray[i]) {
	  fclose(farray[i]);
	}
      }
    }
    fclose(f);
  }

  return 0;

}

