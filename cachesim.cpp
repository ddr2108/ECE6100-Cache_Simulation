#include "cachesim.hpp"
#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <sys/time.h>


//Functions for interacting with cache
void prefetch(uint64_t, cache_stats_t*);
int checkL1(uint64_t, cache_stats_t*, int);
int checkL2(uint64_t, cache_stats_t*, int);
void addToCache(uint64_t, cache_stats_t*, int, int);

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
 * Subroutine that performs prefetching.
 *
 * @address  The target memory address
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

		//Calculate address to prefetch
		address = address + pending_stride*pow(2,b2);
	 	indexForC2 = (address>>b2)&(uint64_t(pow(2,c2-b2-s2)-1));		//Index into cache
		tagForC2 = address>>(c2-s2);										//Tag for cache
		incrementC2 = pow(2,c2-b2-s2);			//Increment between items in same set

		indexOld = indexForC2;
		//Try to find cache slot with oldest item
		for (i=0; i<pow(2,s2); i++){
			//Look for the oldest item in case of removing LRU
			if (validC2[indexForC2 + incrementC2*i] == 1 && (ageC2[indexForC2 + incrementC2*i] < ageC2[indexOld])){
				indexOld = indexForC2 + incrementC2*i;
			}
		}

		//Try to find cache slot with invalid item
		for (i=0; i<pow(2,s2); i++){
			if (validC2[indexForC2 + incrementC2*i] == 0){
				brought = 1;			//if found, mark as a hit

				//Set up timestamp
				ageC2[indexForC2 + incrementC2*i]  = ageC2[indexOld] -1 ;

				//Set up info
				validC2[indexForC2 + incrementC2*i] = 1;
				dirtyC2[indexForC2 + incrementC2*i] = 0;
				//Put tag into cache
				tagC2[indexForC2 + incrementC2*i] = tagForC2;
				break;
			}

		}
			//If all items valid, removed LRU
		if (brought == 0){
			brought = 1;			//if found, mark as a hit

			//Set up timestamp
			ageC2[indexOld]  = ageC2[indexOld] -1 ;

			//Set up info
			validC2[indexOld] = 1;
			dirtyC2[indexOld] = 0;

			//Put tag into cache
			tagC2[indexOld] = tagForC2;
		}
	}


}

/**
 * Subroutine that looks through L1 cache.
 *
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 * @rw The type of event. Either READ or WRITE
 *
 * @return 1 if hit, 0 if not
 */
int checkL1(uint64_t address, cache_stats_t* p_stats, int rw) {
	int i;
	//Flag for hit
	int hit = 0;
	//L1 parameters
	uint64_t indexForC1 = (address>>b1)&(uint64_t(pow(2,c1-b1-s1)-1));		//Index into cache
	uint64_t tagForC1 = address>>(c1-s1);									//Tag for cache
	uint64_t incrementC1 = pow(2,c1-b1-s1);									//Increment between items in same set
	//Struct for time stamps
	struct timeval tv;

	p_stats->L1_accesses++;			//Update stats

	// Check L1 cache - go through each in set
	for (i=0; i<pow(2,s1); i++){
		//Check for valid tag and valididty
		if (tagC1[indexForC1 + incrementC1*i] == tagForC1 && validC1[indexForC1 + incrementC1*i] == 1){
			hit = 1;			//if found, mark as a hit

			//Set up timestamp
			gettimeofday(&tv,NULL);
			ageC1[indexForC1 + incrementC1*i]  = tv.tv_sec*1000000+tv.tv_usec;

			//Change dirty bit based on read or write
			if (rw == WRITE){
				dirtyC1[indexForC1 + incrementC1*i]=1;		//Mark as dirty
			}else{
				dirtyC1[indexForC1 + incrementC1*i]=0;		//Mark as not dirty				
			}
			
			break;		//break if found
		}
	}

	return hit;		//return whether there was a hit
}

/**
 * Subroutine that looks through L2 cache.
 *
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 * @rw The type of event. Either READ or WRITE
 *
 * @return 1 if hit, 0 if not
 */
int checkL2(uint64_t address, cache_stats_t* p_stats, int rw) {
	int i;
	//Flag for hit
	int hit = 0;
	//L2 parameters
	uint64_t indexForC2 = (address>>b2)&(uint64_t(pow(2,c2-b2-s2)-1));		//Index into cache
	uint64_t tagForC2 = address>>(c2-s2);									//Tag for cache
	uint64_t incrementC2 = pow(2,c2-b2-s2);									//Increment between items in same set
	//Struct for time stamps
	struct timeval tv;

	//Update stats
	p_stats->L2_accesses++;			
	if (rw == WRITE){
		p_stats->L1_write_misses++;
	}else{
		p_stats->L1_read_misses++;
	}

	//Check L2 cache - go through each in set
	for (i=0; i<pow(2,s2); i++){
		//Check for valid tag and valididty
		if (tagC2[indexForC2 + incrementC2*i] == tagForC2 && validC2[indexForC2 + incrementC2*i] == 1){
			hit = 2;			//if found, mark as a hit

			//Set up timestamp
			gettimeofday(&tv,NULL);
			ageC2[indexForC2 + incrementC2*i]  = tv.tv_sec*1000000+tv.tv_usec;
				
			//Change dirty bit based on read or write
			if (rw == WRITE){
				dirtyC2[indexForC2 + incrementC2*i]=1;		//Mark as dirty
			}else{
				dirtyC2[indexForC2 + incrementC2*i]=0;		//Mark as not dirty
			}

			break;		//break if found
		}
	}

	return hit;
}

/**
 * Subroutine that looks through L2 cache.
 *
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 * @rw The type of event. Either READ or WRITE
 * @cache The cache. Either L1 or L2
 */
void addToCache(uint64_t address, cache_stats_t* p_stats, int rw, int cache){
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
	int indexOld;			//LRU 
	int brought = 0;			//Flag for bringin into cache
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

	indexOld = indexC;										//Initialize LRU

	//Try to find cache slot with invalid item
	for (i=0; i<pow(2,s); i++){
		//Check if item is invalid
		if (valid[indexC + incrementC*i] == 0){
			brought = 1;			//if found, mark as a hit

			//Set up timestamp
			gettimeofday(&tv,NULL);
			age[indexC + incrementC*i]  = tv.tv_sec*1000000+tv.tv_usec;

			//Set up valid bit
			valid[indexC + incrementC*i] = 1;
			//Change dirty bit based on read or write
			if (rw == WRITE){
				dirty[indexC + incrementC*i] = 1;
			}else{
				dirty[indexC + incrementC*i] = 0;
			}

			//Put tag into cache
			tag[indexC + incrementC*i] = tagC;
				
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

		//Set up timestamp
		gettimeofday(&tv,NULL);
		age[indexOld]  = tv.tv_sec*1000000+tv.tv_usec;

		//Set up valid bit
		valid[indexOld] = 1;
		//Change dirty bit based on read or write
		if (rw == WRITE){
			dirty[indexOld] = 1;
		}else{
			dirty[indexOld] = 0;
		}

		//Put tag into cache
		tag[indexOld] = tagC;
	}

	//Update stats
	if (rw==READ && cache==L2){
		p_stats->L2_read_misses++;		//update stats
	}else if (rw==WRITE && cache==L2){
		p_stats->L2_write_misses++;		//update stats
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
	hit = checkL1(address, p_stats, rw);
	//If there is no hit check L2 cache
	if (hit == 0){
		hit = checkL2(address, p_stats, rw);
	}

	//If there is a miss in the caches
	if (hit != 1){
		//Add to L1 cache
		addToCache(address, p_stats, rw, L1);
		//If miss in L2 cache
		if (hit==0){
			//Add to L2 cache if miss in both cache
			addToCache(address, p_stats, rw, L2);

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
