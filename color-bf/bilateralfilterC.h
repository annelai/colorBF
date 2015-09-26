#ifndef BILATERALFILTER_H
#define BILATERALFILTER_H

#include "fileIO.h"

class BilateralFilterC {
public:
	BilateralFilterC(): r(1), sigma_s(0.33), sigma_r(30), gs_kernel(nullptr), gr_kernel(nullptr) {}
	~BilateralFilterC() {
		if(gs_kernel) {
			delete [] gs_kernel;
		}
		if(gr_kernel) {
			delete [] gr_kernel;
		}
	}
	void initGsGaussian(const int& r);
	void initGrGaussian(const double& sigma_r);
	void directCompute(uchar* I, uchar* O, int width, int height);
    void multiChPBFICMethod(uchar* I, uchar* O, int width, int height, int numSamplePerChannel);
    void multiChPBFICMethod_modified(uchar* I, uchar* O, int width, int height, int numSamplePerChannel);
	void colorSampledPBFICMethod(uchar* I, uchar* O, int* RGB, int width, int height, int numSpace, int numSample);

private:
	int r;
	double sigma_r;
	double sigma_s;
	double* gs_kernel;
	double* gr_kernel;
};




#endif
