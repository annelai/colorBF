#include <iostream>
#include <math.h>
#include "fileIO.h"
#include "timer.h"
#include "bilateralfilterC.h"
using namespace std;

bool verbose = true;

float psnr(uchar* A, uchar* B, int numPixel)
{
	float mse = 0;
	for(int i = 0; i < numPixel; ++i) {
		mse += float((A[i] - B[i])*(A[i] - B[i]));
	}
	mse /= float(numPixel);
	return 20.0*log10(255.0) - 10.0*log10(mse);
}

// Fast bilateral filtering for natural color images
// Usage: <program> <input image> <filter radius> <sigma_r>
int main(int argc, char* argv[])
{
    unsigned seed = (unsigned)time(NULL);
    srand(seed);
	uchar* I;
	uchar* O;
	uchar* O_gold;
    int* RGB;

	int width, height, numPixel, r, numSpace;
	double sigma_s, sigma_r;
	long long tic, toc;

	r = atoi(argv[2]);
	sigma_s = double(r)/3.0;
	sigma_r = atof(argv[3]);


	// Allocate memories & read target images
	sizePGM(width, height, argv[1]);
	numPixel = width*height;
	I = new uchar[numPixel*3];
	O = new uchar[numPixel*3];
	O_gold = new uchar[numPixel*3];
	readPPM(I, argv[1]);

	//verbose && cout << "Input image: " << argv[1] << endl
	//				<< "Size: " << width << "x" << height << endl
	//				<< "filter radius: " << r << endl
	//				<< "sigma_s: " << sigma_s << endl
	//				<< "sigma_r: " << sigma_r << endl;

	BilateralFilterC bf;
	bf.initGsGaussian(r);
    bf.initGrGaussian(sigma_r);

    // Saving filterd image
    string fout = argv[1];
    size_t f1 = fout.find('/') + 1;
    size_t f2 = fout.find('.');
    fout = fout.substr(f1, f2-f1);

	//verbose && cout << "+ BF direct implementation" << endl;
	tic = milliseconds_now();
	memset(O_gold, 0, numPixel*3);
	bf.directCompute(I, O_gold, width, height);
	toc = milliseconds_now();
	//cout << "	Elapsed time: " << float(toc-tic)/1000.f << " sec." << endl;
    //cout << "	Elapsed time: " << float(toc-tic)/(float)CLOCKS_PER_SEC/numPixel*1000000 << " sec." << endl;
	// writePPM(O_gold, width, height, fout + "_cbf_groundtruth_" + argv[3] + ".ppm");

 //    verbose && cout << "+ Multi-channel PBFIC method" << endl;
 //    tic = milliseconds_now();
 //    memset(O, 0, numPixel*3);
 //    bf.multiChPBFICMethod(I, O, width, height, 3);
 //    toc = milliseconds_now();
 //    //cout << "	Elapsed time: " << float(toc-tic)/1000.f << " sec." << endl;
 //    cout << "	Elapsed time: " << float(toc-tic)/(float)CLOCKS_PER_SEC/numPixel*1000000 << " sec." << endl;
	// cout << "	PSNR: " << psnr(O_gold, O, numPixel*3) << " dB." << endl;
 //    writePPM(O, width, height, fout + "_cbf_PBFIC_" + argv[3] + ".ppm");


	// verbose && cout << "+ Color sampled PBFIC method" << endl;
    tic = milliseconds_now();
    sizeColorSpace(numSpace, argv[4]);
    RGB = new int[numSpace*3];
    readColorSpace(RGB, numSpace, argv[4]);
    memset(O, 0, numPixel*3);
    bf.colorSampledPBFICMethod(I, O, RGB, width, height, numSpace, 20);
    toc = milliseconds_now();
    //cout << "	Elapsed time: " << float(toc-tic)/1000.f << " sec." << endl;
    // cout << "	Elapsed time: " << float(toc-tic)/(float)CLOCKS_PER_SEC/numPixel*1000000 << " sec." << endl;
	cout << "	PSNR: " << psnr(O_gold, O, numPixel*3) << " dB." << endl;
    writePPM(O, width, height, fout + "_cbf_KInterpolation_" + argv[3] + ".ppm");


	delete [] I;
	delete [] O;
	delete [] O_gold;
    delete [] RGB;

	return 0;

}
