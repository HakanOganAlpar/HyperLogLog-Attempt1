#pragma once

#ifndef HYPERLOGLOGH
#define HYPERLOGLOGH


#define WINDOWSIZE 4;
#define HASHSEED 7;


#include <string>
#include <iostream>
#include <vector>
using namespace std;

class HyperLogLog
{
public:
	HyperLogLog(unsigned int bitsTaken1);

	void Add (char param [], int length);

	void Print();

	void EstimationEQ();

private:
	unsigned int  numBuckets;
	int hashSeed = HASHSEED;
	int** buckets;
	int windowSize = WINDOWSIZE;
	int bitsTaken;
	unsigned int Encode(char * p);
	int cardinality;

};













#endif
