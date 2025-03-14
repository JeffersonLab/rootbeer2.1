#include "particleDEF.h"
#define PI  		3.14159
#define RAD2DEG  	(180.0/PI)

typedef struct count_part_t {
        float   pid;
        float   x;
        float   y;
	float   z;
	float   E;
	float   px;
	float   py;
	float   pz;
	float   q;
	int     trkid;
	int     sec;
	float   beta;
	int     sc_stat;
	float   sc_time;
	float   sc_len;
	int     st_stat;
	float   st_time;
	float   st_len;
	float   mass;
	int     mass_ref;
	float   betam;
	float   epho;
	float   tpho;
	int     g_stat;
	int     ngrf;
} count_part_t;
