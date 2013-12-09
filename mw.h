#ifndef MULTI_WIGGLE
#define MULTI_WIGGLE

#define VALUE float

typedef unsigned long ULONG;
typedef unsigned short USHORT;
typedef unsigned int UINT;
typedef unsigned char UCHAR;

#define MAX_TRACKS 50
#define FORMAT_MWIG "wig"
#define FORMAT_MBED "bed"
#define VERSION_MAJOR 0
#define VERSION_MINOR 2

#define META_KEY_LEN 32
#define META_VALUE_LEN 255
#define REGION_NAME_LEN 32

// meta key from the header
typedef struct {
  char k[META_KEY_LEN];
  char v[META_VALUE_LEN];
} META;

// track record
typedef struct {
  USHORT id;
  VALUE min;
  VALUE max;
  UCHAR ori;
  char name[32];
  char desc[255];
} TRACK; 

// region record
typedef struct {
  char name[REGION_NAME_LEN];
  ULONG offset;
  ULONG size;
} REGION;

typedef struct {
  TRACK t;
  VALUE *v;
} RESULT;


#define META_NUM 10 // fields in the header

typedef struct {
  char format[5]; // MWIG or MBED ? 
  USHORT v_major; // sw and data should macth to be compatible
  USHORT v_minor; // minor changes not affecting file format
  USHORT tracks; // number of tracks in the file
  USHORT regions; // number of regions in the file
  ULONG offset; // where the real data start
  UINT taxon; // species id, e.g 284812 for s.pombe
  char assembly[64]; // assembly id in ENA, e.g GCA_000002945.2
  VALUE min; //min value across all tracks
  VALUE max; //max value across all tracks
  char desc[255]; //description of the file, e.g study GSE24360
} HEADER;


typedef struct T {
  TRACK data;
  char name[32];
  struct T* prev;
  struct T* next;
} TRECORD;


typedef struct R {
  REGION data;
  struct R* prev;
  struct R* next;
} RRECORD;


int mw_create(char *fname, char **argv, int argc);
META *mw_stats(char *fname);
TRACK *mw_tracks (char *fname);
REGION *mw_regions(char *fname);
RESULT *mw_fetch(char *fname, char *region, char *tracks, int winsize, int *tcount);
int mw_dump(char *fname, char *region, char *tracks);
 
#endif
