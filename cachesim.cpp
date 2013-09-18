#include "cachesim.hpp"
#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <sys/time.h>


//Functions for interacting with cache
void prefetch(uint64_t, cache_stats_t*);
int checkCache(uint64_t, cache_stats_t*, int, int);
void addToCache(uint64_t, cache_stats_t*, int, int, int);

//Cache pieces
uint64_t* tagC1;		//Store Tags 
long long* ageC1;			//Store age via timestamps
int* validC1;			//Information bits
int* dirtyC1;
uint64_t* tagC2;		//Store Tags 
long long* ageC2;			//Store age via timestamps
int* validC2;			//Information bits
int* dirtyC2;
//Cache stats
uint64_t c1;
uint64_t b1;
uint64_t s1;
uint64_t c2;
uint64_t b2;
uint64_t s2;
uint32_t k;
//Variables for prefetching
uint64_t prev_miss = -1;
uint64_t pending_stride = -1;

//Cache action
#define READ  'r'
#define WRITE 'w'
//Which cache
#define L1  1
#define L2  2
//Prefectch
#define NO_PREFETCH  0
#define PREFETCH     1
/**
 * Subroutine for initializing the cache. You many add and initialize any global or heap
 * variables as needed.
 *
 * @c1 Total size of L1 in bytes is 2^C1
 * @b1 Size of each block in L1 in bytes is 2^B1
 * @s1 Number of blocks per set in L1 is 2^S1
 * @c2 Total size of L2 in bytes is 2^C2
 * @b2 Size of each block in L2 in bytes is 2^B2
 * @s2 Number of blocks per set in L2 is 2^S2
 * @k Prefetch K subsequent blocks
 */
void setup_cache(uint64_t in_c1, uint64_t in_b1, uint64_t in_s1, uint64_t in_c2, uint64_t in_b2, uint64_t in_s2, uint32_t in_k) {
	int i;
	
	//Set stats in global 
	c1 = in_c1;
	b1 = in_b1;
	s1 = in_s1;
	c2 = in_c2;
	b2 = in_b2;
	s2 = in_s2;
	k = in_k;

	//Allocate memory for L1 cache
	tagC1 = (uint64_t*) malloc((int)pow(2,(c1-b1))*sizeof(uint64_t));
	ageC1 = (long long*) malloc((int)pow(2,(c1-b1))*sizeof(long long));
	validC1 = (int*) malloc((int)pow(2,(c1-b1))*sizeof(int));
	dirtyC1 = (int*) malloc((int)pow(2,(c1-b1))*sizeof(int));
	//Allocate memory for L2 cache
	tagC2 = (uint64_t*) malloc((int)pow(2,(c2-b2))*sizeof(uint64_t));
	ageC2 = (long long*) malloc((int)pow(2,(c2-b2))*sizeof(long long));
	validC2 = (int*) malloc((int)pow(2,(c2-b2))*sizeof(int));
	dirtyC2 = (int*) malloc((int)pow(2,(c2-b2))*sizeof(int));


	//Initialize valid bits to 0
	for (i = 0; i<(int)pow(2,(c1-b1)); i++){
		validC1[i] = 0;
	}
	for (i = 0; i<(int)pow(2,(c2-b2)); i++){
		validC2[i] = 0;
	}

}

/**
 * Subroutine that looks through caches.
 *
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 * @rw The type of event. Either READ or WRITE
 * @cache The cache. Either L1 or L2
 *
 * @return 1 if hit, 0 if not
 */
int checkCache(uint64_t address, cache_stats_t* p_stats, int rw, int cache) {
	int i;
	//Flag for hit
	int hit = 0;
	//Cache parameters
	int c,b,s;
	uint64_t* tag;
	long long* age;	
	int* valid, *dirty;
	//Address parameters
	uint64_t indexC;					//Index into cache
	uint64_t tagC;						//Tag for cache
	uint64_t incrementC;				//Increment between items in same set
	//Struct for time stamps
	struct timeval tv;

	//Adjust cache parameters based on the cache
	if (cache==L1){
		c = c1;
		b = b1;
		s = s1;
		valid = validC1;
		dirty = dirtyC1;
		age = ageC1;
		tag = tagC1;
	} else{
		c = c2;
		b = b2;
		s = s2;
		valid = validC2;
		dirty = dirtyC2;
		age = ageC2;
		tag = tagC2;
	}

	//Determine address parameters
	indexC = (address>>b)&(uint64_t(pow(2,c-b-s)-1));		//Index into cache
	tagC = address>>(c-s);									//Tag for cache
	incrementC = pow(2,c-b-s);								//Increment between items in same set

	// Check cache - go through each in set
	for (i=0; i<pow(2,s); i++){
		//Check for valid tag and valididty
		if (tag[indexC + incrementC*i] == tagC && valid[indexC + incrementC*i] == 1){
			//Adjust value of hit based on which cache missed
			if (cache==L1){
				hit = 1;
			}else if (cache==L2){
				hit = 2;
			}
			//Set up timestamp
			gettimeofday(&tv,NULL);
			age[indexC + incrementC*i]  = tv.tv_sec*1000000+tv.tv_usec;
			//Change dirty bit based on read or write
			if (rw == WRITE){
				dirty[indexC + incrementC*i]=1;		//Mark as dirty
			}else{
				dirty[indexC + incrementC*i]=0;		//Mark as not dirty				
			}	
			break;		//break if found
		}
	}

	//Update access stats
	if (cache==L1){
		p_stats->L1_accesses++;
	}else if (cache==L2){
		p_stats->L2_accesses++;
	}
	//Update misses stats
	if (cache == L1 && rw==READ && hit == 0){
		p_stats->L1_read_misses++;
	}else if (cache == L2 && rw==READ && hit == 0){
		p_stats->L2_read_misses++;
	} else if (cache == L1 && rw==WRITE && hit == 0){
		p_stats->L1_write_misses++;
	}else if (cache == L2 && rw==WRITE && hit == 0){
		p_stats->L2_write_misses++;
	}

	return hit;		//return whether there was a hit
}

/**
 * Subroutine that looks through L2 cache.
 *
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 * @rw The type of event. Either READ or WRITE
 * @cache The cache. Either L1 or L2
 */
void addToCache(uint64_t address, cache_stats_t* p_stats, int rw, int cache, int prefetch){
	int i;
	//Cache parameters
	int c,b,s;
	uint64_t* tag;
	long long* age;	
	int* valid, *dirty;
	//Address parameters
	uint64_t indexC;					//Index into cache
	uint64_t tagC;						//Tag for cache
	uint64_t incrementC;				//Increment between items in same set
	//Cache Access
	int indexOld;				//LRU 
	int indexPrefetch;
	int brought = 0;			//Flag for bringin into cache
	int indexAdded;
	//Struct for time stamps
	struct timeval tv;

	//Adjust cache parameters based on the cache
	if (cache==L1){
		c = c1;
		b = b1;
		s = s1;
		valid = validC1;
		dirty = dirtyC1;
		age = ageC1;
		tag = tagC1;
	} else{
		c = c2;
		b = b2;
		s = s2;
		valid = validC2;
		dirty = dirtyC2;
		age = ageC2;
		tag = tagC2;
	}

	//Determine address parameters
	indexC = (address>>b)&(uint64_t(pow(2,c-b-s)-1));		//Index into cache
	tagC = address>>(c-s);									//Tag for cache
	incrementC = pow(2,c-b-s);								//Increment between items in same set

	//Find oldest item if prefetching
	if (prefetch == PREFETCH){
		indexPrefetch = indexC;
		//Try to find cache slot with oldest item
		for (i=0; i<pow(2,s); i++){
			//Look for the oldest item for timing purposes
			if (valid[indexC + incrementC*i] == 1 && (age[indexC + incrementC*i] < age[indexPrefetch])){
				indexPrefetch = indexC + incrementC*i;
			}
		}
	}

	indexOld = indexC;				//Initialize LRU

	//Try to find cache slot with invalid item
	for (i=0; i<pow(2,s); i++){
		//Check if item is invalid
		if (valid[indexC + incrementC*i] == 0){
			brought = 1;			//if found, mark as a hit
			//Index added to
			indexAdded = indexC + incrementC*i;		
			break;
		}

		//Look for the oldest item in case of removing LRU
		if (age[indexC + incrementC*i] < age[indexOld]){
			indexOld = indexC + incrementC*i;
		}
	}
	
	//If all items valid, removed LRU
	if (brought == 0){
		brought = 1;					//set flag that brough in
		//Check if need to do writeback
		if (dirty[indexOld] == 1){
			p_stats->write_backs++;
		}
		//Index to add to
		indexAdded = indexOld;
	}

	//Set up data for added item
	//Set up timestamp
	if (prefetch == PREFETCH){
		//If prefetch set timestamp to older than oldest
		age[indexAdded] = age[indexPrefetch] - 1;
	}else{
		//Set up timestamp
		gettimeofday(&tv,NULL);
		age[indexAdded]  = tv.tv_sec*1000000+tv.tv_usec;
	}
	//Set up valid bit
	valid[indexAdded] = 1;
	//Change dirty bit based on read or write
	if (rw == WRITE){
		dirty[indexAdded] = 1;
	}else{
		dirty[indexAdded] = 0;
	}
	//Put tag into cache
	tag[indexAdded] = tagC;
}


/**
 * Subroutine that performs prefetching.
 *
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 */
void prefetch(uint64_t address, cache_stats_t* p_stats){
	int i,j = 0;
	int64_t d = 0;		//stride between misses
	int brought = 0;	//flag for bringin in
	//L2 parameters
	uint64_t indexForC2 = (address>>b2)&(uint64_t(pow(2,c2-b2-s2)-1));		//Index into cache
	uint64_t tagForC2 = address>>(c2-s2);										//Tag for cache
	uint64_t incrementC2 = pow(2,c2-b2-s2);			//Increment between items in same set

	//For storing LRU
	int indexOld = 0;

	//Store initial miss if never missed before
	if (prev_miss == -1){
		prev_miss = (address>>b2);
		return;
	}

	//Calculate stride
	d = (address>>b2) - prev_miss;

	prev_miss = address>>b2;		//store new miss
	if (d != pending_stride){
		pending_stride = d;
		return;
	}
	pending_stride = d;				//set stride

	//Perform prefetching
	for (j=1; j<=k; j++){
		p_stats->prefetched_blocks++;	//update stats
				address = address + pending_stride*pow(2,b2);

		//Add to L2 cache for prefetching
		addToCache(address, p_stats, READ, L2, PREFETCH);
	}


}

/**
 * Subroutine that simulates the cache one trace event at a time.
 *
 * @rw The type of event. Either READ or WRITE
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 */
void cache_access(char rw, uint64_t address, cache_stats_t* p_stats) {
	int hit = 0;			//flag for hits

	//Increment number of accesses
	p_stats->accesses++;		
	//Check if read or write was done
	if (rw==READ){
		p_stats->reads++;
	}else{
		p_stats->writes++;
	}

	//Check if address in L1 cache
	hit = checkCache(address, p_stats, rw, L1);
	//If there is no hit check L2 cache
	if (hit == 0){
		hit = checkCache(address, p_stats, rw, L2);
	}

	//If there is a miss in the caches
	if (hit != 1){
		//Add to L1 cache
		addToCache(address, p_stats, rw, L1, NO_PREFETCH);
		//If miss in L2 cache
		if (hit==0){
			//Add to L2 cache if miss in both cache
			addToCache(address, p_stats, rw, L2, NO_PREFETCH);

			//Prefetch anything else
			prefetch(address, p_stats);
		}
	}
}


/**
 * Subroutine for cleaning up any outstanding memory operations and calculating overall statistics
 * such as miss rate or average access time.
 *
 * @p_stats Pointer to the statistics structure
 */
void complete_cache(cache_stats_t *p_stats) {
	double HT1, HT2 = 0;			//Variables for AAT
	double MP1 = 0;
	double MP2 = 500.0;
	double MR1, MR2 = 0;

	//Free memory associated with L1 cache
	free(tagC1);
	free(ageC1);
	free(validC1);
	free(dirtyC1);
	//Free memory associated with L2 cache
	free(tagC2);
	free(ageC2);
	free(validC2);
	free(dirtyC2);

	//Calculate AAT
	MR1 = (double)(p_stats->L1_read_misses + p_stats->L1_write_misses)/p_stats->L1_accesses;
	MR2 = (double)(p_stats->L2_read_misses + p_stats->L2_write_misses)/p_stats->L2_accesses;
	HT1 = 2 + 0.2*s1;
	HT2 = 4 + 0.4*s2;
	MP1 = HT2 + MR2*MP2;
	p_stats->avg_access_time = HT1 + MR1 * MP1;
}
