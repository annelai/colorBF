#include "poissondisk.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

void updatequeue(vec3<int>* p, int& queueSize, int* queue, vec3<int>* center, int& sampleSize, float r, bool all) {
    // update valid queue.
    if (all) {
        for (int k = 0; k < sampleSize; ++k)
        {
            for (int q = 0; q < queueSize; ++q) {
                float dist = (center[k].x-p[queue[q]].x)*(center[k].x-p[queue[q]].x) \
                       + (center[k].y-p[queue[q]].y)*(center[k].y-p[queue[q]].y) \
                       + (center[k].z-p[queue[q]].z)*(center[k].z-p[queue[q]].z);
                if (dist < r) {
                    queue[q] = queue[queueSize-1];
                    --q;
                    --queueSize;
                }
            }
        }
    }
    else {
        for (int q = 0; q < queueSize; ++q) {
            float dist = (center[sampleSize-1].x-p[queue[q]].x)*(center[sampleSize-1].x-p[queue[q]].x) \
                   + (center[sampleSize-1].y-p[queue[q]].y)*(center[sampleSize-1].y-p[queue[q]].y) \
                   + (center[sampleSize-1].z-p[queue[q]].z)*(center[sampleSize-1].z-p[queue[q]].z);
            if (dist < r) {
                queue[q] = queue[queueSize-1];
                --q;
                --queueSize;
            }
        }
    }
}

// Bridson’s Algorithm
void PoissonDisk(int* Ir, int* Ig, int* Ib, vec3<int>* center, int numSpace, int numSample, float r)
{
    //cout << "numSpace: " << numSpace << endl;
    int sampleSize = 0;
    int queueSize = numSpace;
    int* queue = new int[numSpace];
    vec3<int>* p = new vec3<int>[numSpace];
    for (int i = 0; i < numSpace; ++i) {
        p[i].x = Ir[i];
        p[i].y = Ig[i];
        p[i].z = Ib[i];
        queue[i] = i;
    }
    // Initialize a sample.
    int s = queue[rand() % queueSize];
    center[0].x = p[s].x;
    center[0].y = p[s].y;
    center[0].z = p[s].z;
    ++sampleSize;
    queue[s] = queue[queueSize-1];
    --queueSize;

    while(sampleSize < numSample && queueSize) {
        updatequeue(p, queueSize, queue, center, sampleSize, r, 0);
        // randomly pick a sample from the pool.
        int s = queue[rand() % queueSize];
        center[sampleSize].x = p[queue[s]].x;
        center[sampleSize].y = p[queue[s]].y;
        center[sampleSize].z = p[queue[s]].z;
        ++sampleSize;
        queue[s] = queue[queueSize-1];
        --queueSize;
    }
    //cout << queueSize << " left in queue." << endl;
}

void AutoPoissonDisk(int* Ir, int* Ig, int* Ib, vec3<int>* center, int numSpace, int numSample) {

    //cout << "numSpace: " << numSpace << endl;
    float r = 1500;
    float rate = 0.9;

    int sampleSize = 0;
    int queueSize = numSpace;
    int* queue = new int[numSpace];
    vec3<int>* p = new vec3<int>[numSpace];
    for (int i = 0; i < numSpace; ++i) {
        p[i].x = Ir[i];
        p[i].y = Ig[i];
        p[i].z = Ib[i];
        queue[i] = i;
    }
    // Initialize a sample.
    int s = queue[rand() % queueSize];
    center[0].x = p[s].x;
    center[0].y = p[s].y;
    center[0].z = p[s].z;
    ++sampleSize;
    queue[s] = queue[queueSize-1];
    --queueSize;

    while(sampleSize < numSample && r > 50) {
        updatequeue(p, queueSize, queue, center, sampleSize, r, 0);
        if (queueSize == 0) {
            r *= rate;
            queueSize = numSpace;
            for (int i = 0; i < numSpace; ++i) queue[i] = i;
            updatequeue(p, queueSize, queue, center, sampleSize, r, 1);
        }
        else {
            // randomly pick a sample from the pool.
            s = queue[rand() % queueSize];
            center[sampleSize].x = p[queue[s]].x;
            center[sampleSize].y = p[queue[s]].y;
            center[sampleSize].z = p[queue[s]].z;
            ++sampleSize;
            queue[s] = queue[queueSize-1];
            --queueSize;

        }
    }
    //cout << queueSize << " left in queue." << endl;
}


