#include "bilateralfilterC.h"
#include "matrix.h"
#include "fileIO.h"
#include "poissondisk.h"
#include <cmath>

void box_filter(double* I, double* Out, int width, int height, int r) {
    for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
        float acc = 0;
        for (int r1 = -r; r1 < r; ++r1) {
        for (int r2 = -r; r2 < r; ++r2) {
            if ( (y+r1) >= 0 && (y+r1) < height && (x+r2) >= 0 && (x+r2) < width ) {
                acc += float(I[(y+r1)*width+x+r2]);
            }
        }
        }
        Out[y*width+x] = acc;
    }
    }
}

void box_filter_1(double* I, double* O, int width, int height, int r) {
    double* sum = new double[width*height];
    double* weight = new double[width*height];
    int min_i, min_j;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {

            min_i = max(i-1, 0);
            min_j = max(j-1, 0);
            sum[i*width+j] = I[i*width+j] + sum[min_i*width+j] - sum[min_i*width+min_j] + sum[i*width+min_j];
            weight[i*width+j] = 1 + weight[min_i*width+j] - weight[min_i*width+min_j] + weight[i*width+min_j];
        }
    }
    int max_y, min_y, max_x, min_x;
    double w;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {

            max_y = min(y+r, height-1);
            min_y = max(y-r-1, 0);
            max_x = min(x+r, width-1);
            min_x = max(x-r-1, 0);

            O[y*width+x] = sum[max_y*width+max_x] - sum[min_y*width+max_x] - sum[max_y*width+min_x] + sum[min_y*width+min_x];
            w = weight[max_y*width+max_x] - weight[min_y*width+max_x] - weight[max_y*width+min_x] + weight[min_y*width+min_x];

            O[y*width+x] /= w;
        }
    }

}

void box_filter_2(double* I, double* O, int width, int height, int r) {
    int len = 2*r + 1;
    double* hist_acc = new double[len];
    double* hist_wacc = new double[len];

    int x_begin, x_end, y_begin, y_end;
    double acc, wacc;
    for (int y = 0; y < height; ++y) {
        // x = 0
        y_begin = max(0, y-r);
        x_begin = 0;
        y_end = min(height-1, y+r);
        x_end = min(width-1, r);
        // Reset the histogram
        for(int s = 0; s < len; ++s) {
            hist_acc[s] = 0;
            hist_wacc[s] = 0;
        }

        // Initialize the histogram
        for(int sx = x_begin; sx <= x_end; ++sx) {
            for(int sy = y_begin; sy <= y_end; ++sy) {
                hist_acc[sx+r] += I[sy*width+sx];
                ++hist_wacc[sx+r];
            }
        }
        // Aggregate the weighted sum
        acc = wacc = 0;
        for(int s = 0; s < len; ++s) {
            if(hist_wacc[s] > 0) {
                acc += hist_acc[s];
                wacc += hist_wacc[s];
            }
        }
        O[y*width] = acc/wacc;
        // x > 0, the histogram can be shared
        for(int x = 1; x < width; ++x) {
            x_begin = max(0, x-r);
            x_end = min(width-1, x+r);

            // Shift the histogram
            for(int s = 0; s < len-1; ++s) {
                --hist_acc[s];
                hist_wacc[s] = hist_wacc[s+1];
            }
            // the right most position
            hist_acc[len-1] = 0;
            hist_wacc[len-1] = 0;
            // Add the new column
            for(int sy = y_begin; sy <= y_end; ++sy) {
                hist_acc[len-1] += I[sy*width+x_end];
                ++hist_wacc[len-1];
            }
            // Aggregate the weighted sum
            acc = wacc = 0;
            for(int s = 0; s < len; ++s) {
                if(hist_wacc[s] > 0) {
                    acc += hist_acc[s];
                    wacc += hist_wacc[s];
                }
            }
            O[y*width+x] = acc/wacc;
        } //end x
    } // end y
}

void boxfilter(double* I, double* Out, int width, int height, int r, double* tmp)
{
    // cumulative sum over Y axis
    for(int x = 0; x < width; ++x) {
        // y = 0
        tmp[x] = I[x];
        // y > 0
        for(int y = 1; y < height; ++y) {
            tmp[y*width+x] = tmp[(y-1)*width+x] + I[y*width+x];
        }
    }

    // difference over Y axis
    // imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
    for(int y = 0; y <= r; ++y) {
    for(int x = 0; x < width; ++x) {
        Out[y*width+x] = tmp[(y+r)*width+x];
    }
    }
    // imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
    for(int y = r+1; y < height-r; ++y) {
    for(int x = 0; x < width; ++x) {
        Out[y*width+x] = tmp[(y+r)*width+x] - tmp[(y-r-1)*width+x];
    }
    }
    // imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);
    for(int y = height-r; y < height; ++y) {
    for(int x = 0; x < width; ++x) {
        Out[y*width+x] = tmp[(height-1)*width+x] - tmp[(y-r-1)*width+x];
    }
    }

    // cumulative sum over X axis
    for(int y = 0; y < height; ++y) {
        // x = 0
        tmp[y*width] = Out[y*width];
        // x > 0
        for(int x = 1; x < width; ++x) {
            tmp[y*width+x] = tmp[y*width+x-1] + Out[y*width+x];
        }
    }
    // difference over X axis
    // imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
    for(int x = 0; x <= r; ++x) {
    for(int y = 0; y < height; ++y) {
        Out[y*width+x] = tmp[y*width+x+r];
    }
    }
    // imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
    for(int x = r+1; x < width-r; ++x) {
    for(int y = 0; y < height; ++y) {
        Out[y*width+x] = tmp[y*width+x+r] - tmp[y*width+x-r-1];
    }
    }
    // imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);
    for(int x = width-r; x < width; ++x) {
    for(int y = 0; y < height; ++y) {
        Out[y*width+x] = tmp[y*width+width-1] - tmp[y*width+x-r-1];
    }
    }
}

void gaussian1D_recursive(double* I, double* Out, int width, int height, double sigma, double* tmp1, double* tmp2, double* tmp) {

	double a0, a1, a2, a3, b1, b2;
	double common_term1 = exp(-1.695f/sigma);
	b2 = exp(-3.39f/sigma);
	a0 = (1.f - common_term1)*(1.f - common_term1)/(1.f + 3.39f*common_term1/sigma - b2);
	b1 = -2.f*common_term1;
	a1 = (1.695f/sigma - 1.f)*common_term1*a0;
	a2 = (1.695f/sigma + 1.f)*common_term1*a0;
	a3 = -a0*b2;
	// Filtering along horizontal direction
	for(int y = 0; y < height; ++y) {
		// Initial condition: y0 = a0x0
		tmp1[y*width] = a0*double(I[y*width]);
		tmp1[y*width+1] = a0*double(I[y*width+1]) + a1*double(I[y*width]) - b1*tmp1[y*width];
		for(int x = 2; x < width ; ++x) {
			tmp1[y*width+x] = a0*double(I[y*width+x]) + a1*double(I[y*width+x-1]) - b1*tmp1[y*width+x-1] - b2*tmp1[y*width+x-2];
		}
		// Anti-causal filtering
		tmp2[y*width+width-1] = a2*double(I[y*width+width-1]);
		tmp2[y*width+width-2] = a2*double(I[y*width+width-1]) - b1*tmp2[y*width+width-1];
		for(int x = width-3; x >= 0; --x) {
			tmp2[y*width+x] = a2*double(I[y*width+x+1]) + a3*double(I[y*width+x+2]) - b1*tmp2[y*width+x+1] - b2*tmp2[y*width+x+2];
		}
	}
	for(int i = 0; i < height*width; ++i) {
		tmp[i] = tmp1[i] + tmp2[i];
	}

	// Filtering along vertical direction
	for(int x = 0; x < width; ++x) {
		// Initial condition: y0 = a0x0
		tmp1[x] = a0*tmp[x];
		tmp1[width+x] = a0*tmp[width+x] + a1*tmp[x] - b1*tmp1[x];
		for(int y = 2; y < height ; ++y) {
			tmp1[y*width+x] = a0*tmp[y*width+x] + a1*tmp[(y-1)*width+x] - b1*tmp1[(y-1)*width+x] - b2*tmp1[(y-2)*width+x];
		}
		// Anti-causal filtering
		tmp2[(height-1)*width+x] = a2*tmp[(height-1)*width+x];
		tmp2[(height-2)*width+x] = a2*tmp[(height-1)*width+x] - b1*tmp2[(height-1)*width+x];
		for(int y = height-3; y >= 0; --y) {
			tmp2[y*width+x] = a2*tmp[(y+1)*width+x] + a3*tmp[(y+2)*width+x] - b1*tmp2[(y+1)*width+x] - b2*tmp2[(y+2)*width+x];
		}
	}
	for(int i = 0; i < height*width; ++i) {
		Out[i] = tmp1[i] + tmp2[i];
	}
}

/*
 * @param r Filter radius; spatial filter strength sigma_s will be double(r)/3
 * @param sigma_r Range filter strength
*/
void BilateralFilterC::initGsGaussian(const int& r1)
{
	int len = 2*r1 + 1;
	r = r1;
	sigma_s = double(r1)/3.0;

	// Precompute the spatial Gaussian weights
	gs_kernel = new double[len];
	for(int x = -r; x <= r; ++x) {
		gs_kernel[x+r] = exp(-0.5f*double(x*x)/sigma_s/sigma_s);
	}
}

void BilateralFilterC::initGrGaussian(const double& sigma_r1)
{
	sigma_r = sigma_r1;
	// Precompute the range Gaussian weights
	gr_kernel = new double[256];
	for(int i = 0; i < 256; ++i) {
		gr_kernel[i] = exp(-0.5f*double(i*i)/sigma_r/sigma_r);
	}
}

/** Direct implementation of bilateral filter
 * @param I Input image
 * @param O Output image
 * @param width Image width
 * @param height Image height
*/
void BilateralFilterC::directCompute(uchar* I, uchar* O, int width, int height)
{
	int numPixel = width*height;
	uchar* Ir = new uchar[numPixel];
	uchar* Ig = new uchar[numPixel];
	uchar* Ib = new uchar[numPixel];

	for(int i = 0; i < numPixel; ++i) {
		Ir[i] = I[3*i+0];
		Ig[i] = I[3*i+1];
		Ib[i] = I[3*i+2];
	}

	int x_begin, x_end, y_begin, y_end;
	double acc_r, acc_g, acc_b, wacc;
	for(int y = 0; y < height; ++y) {
	for(int x = 0; x < width; ++x) {
		y_begin = max(0, y-r);
		x_begin = max(0, x-r);
		y_end = min(height-1, y+r);
		x_end = min(width-1, x+r);
		acc_r = acc_g = acc_b = wacc = 0;
		for(int sy = y_begin; sy <= y_end; ++sy) {
		for(int sx = x_begin; sx <= x_end; ++sx) {
			double gr = gr_kernel[abs(Ir[sy*width+sx]-Ir[y*width+x])]*gr_kernel[abs(Ig[sy*width+sx]-Ig[y*width+x])]*gr_kernel[abs(Ib[sy*width+sx]-Ib[y*width+x])];
            double gs = 1.0;
			// double gs = gs_kernel[sy-y+r]*gs_kernel[sx-x+r];
			acc_r += gs*gr*Ir[sy*width+sx];
			acc_g += gs*gr*Ig[sy*width+sx];
			acc_b += gs*gr*Ib[sy*width+sx];
			wacc += gs*gr;
		}
		}
		O[3*(y*width+x)  ] = checkPixelRange(acc_r/wacc);
		O[3*(y*width+x)+1] = checkPixelRange(acc_g/wacc);
		O[3*(y*width+x)+2] = checkPixelRange(acc_b/wacc);
	}
	}

	delete [] Ir;
	delete [] Ig;
	delete [] Ib;
}


/** Implementation of Qingxiong Yang's work.
 *  "Constant Time Median and Bilateral Filtering", appeared in Springer 2014.
 *  The method is an expansion of the idea based on PBFIC method, which reduces computation by only conducting a few linear filtering and interpolates the others.
 * @param I Input image
 * @param O Output image
 * @param width Image width
 * @param height Image height
 * @param num The number of PBFIC to be used in each channel
*/

void BilateralFilterC::multiChPBFICMethod(uchar* I, uchar* O, int width, int height, int num)
{
    int numPixel = height*width;
    uchar* Ir = new uchar[numPixel];
    uchar* Ig = new uchar[numPixel];
    uchar* Ib = new uchar[numPixel];
    for(int i = 0; i < numPixel; ++i) {
        Ir[i] = I[3*i+0];
        Ig[i] = I[3*i+1];
        Ib[i] = I[3*i+2];
    }

    // Initialize uniform PBFIC
    int* value = new int[num];
    for(int i = 0; i < num; ++i) {
        value[i] = (255*i)/(num-1);
    }
    // Filtered components
    double** Jr = new double*[num*num*num];
    double** Jg = new double*[num*num*num];
    double** Jb = new double*[num*num*num];
    double* Jw = new double[numPixel];
    // Some tmp buffers
    double* tmp  = new double[numPixel];
    double* tmp_r  = new double[numPixel];
    double* tmp_g  = new double[numPixel];
    double* tmp_b  = new double[numPixel];
    double* tmp1 = new double[numPixel];
    double* tmp2 = new double[numPixel];
    double* tmp3 = new double[numPixel];

    // Compute bilateral filtering for each color exemplar
    for (int i = 0; i < num; ++i) {
    for (int j = 0; j < num; ++j) {
    for (int k = 0; k < num; ++k) {
        Jr[num*num*i+num*j+k] = new double[numPixel];
        Jg[num*num*i+num*j+k] = new double[numPixel];
        Jb[num*num*i+num*j+k] = new double[numPixel];
        for(int pos = 0; pos < numPixel; ++pos) {
            tmp[pos] = gr_kernel[abs(value[i]-Ir[pos])]*gr_kernel[abs(value[j]-Ig[pos])]*gr_kernel[abs(value[k]-Ib[pos])];
        }
        gaussian1D_recursive(tmp, Jw, width, height, sigma_s, tmp1, tmp2, tmp3);
        // boxfilter(tmp, Jw, width, height, r, tmp1);
        for(int pos = 0; pos < numPixel; ++pos) {
            tmp_r[pos] = tmp[pos]*double(Ir[pos]);
            tmp_g[pos] = tmp[pos]*double(Ig[pos]);
            tmp_b[pos] = tmp[pos]*double(Ib[pos]);
        }
        gaussian1D_recursive(tmp_r, Jr[num*num*i+num*j+k], width, height, sigma_s, tmp1, tmp2, tmp3);
        gaussian1D_recursive(tmp_g, Jg[num*num*i+num*j+k], width, height, sigma_s, tmp1, tmp2, tmp3);
        gaussian1D_recursive(tmp_b, Jb[num*num*i+num*j+k], width, height, sigma_s, tmp1, tmp2, tmp3);
        // boxfilter(tmp_r, Jr[num*num*i+num*j+k], width, height, r, tmp1);
        // boxfilter(tmp_g, Jg[num*num*i+num*j+k], width, height, r, tmp1);
        // boxfilter(tmp_b, Jb[num*num*i+num*j+k], width, height, r, tmp1);
        for(int pos = 0; pos < numPixel; ++pos) {
            Jr[num*num*i+num*j+k][pos] /= Jw[pos];
            Jg[num*num*i+num*j+k][pos] /= Jw[pos];
            Jb[num*num*i+num*j+k][pos] /= Jw[pos];
        }
    }
    }
    }
    // Output selection
    bool isDone;
    for(int pos = 0; pos < numPixel; ++pos) {
        isDone = 0;
        for (int i = 0; i < num; ++i) {
            for (int j = 0; j < num; ++j) {
                for (int k = 0; k < num; ++k) {
                    if (Ir[pos] == value[i] && Ig[pos] == value[j] && Ib[pos] == value[k] ) {
                        O[3*pos  ] = checkPixelRange(Jr[num*num*i+num*j+k][pos]);
                        O[3*pos+1] = checkPixelRange(Jg[num*num*i+num*j+k][pos]);
                        O[3*pos+2] = checkPixelRange(Jb[num*num*i+num*j+k][pos]);
                        isDone = 1;
                        break;
                    }
                }
                if(isDone) break;
            }
            if(isDone) break;
        }
        if (isDone) continue;

        int idx_r=0, idx_g=0, idx_b=0;
        double a1=0, a2=0, a3=0, b1=0, b2=0, b3=0;
        for (int i = 0; i < num; ++i) {
            double L = value[i+1]-value[i];
            if(Ir[pos] >= value[i] && Ir[pos] <= value[i+1]) {
                idx_r = i;
                a1 = (value[i+1]-Ir[pos])/L;
                b1 = (Ir[pos] - value[i])/L;
                break;
            }
        }
        for (int i = 0; i < num; ++i) {
            double L = value[i+1]-value[i];
            if(Ig[pos] >= value[i] && Ig[pos] <= value[i+1]) {
                idx_g = i;
                a2 = (value[i+1]-Ig[pos])/L;
                b2 = (Ig[pos] - value[i])/L;
                break;
            }
        }
        for (int i = 0; i < num; ++i) {
            double L = value[i+1]-value[i];
            if(Ib[pos] >= value[i] && Ib[pos] <= value[i+1]) {
                idx_b = i;
                a3 = (value[i+1]-Ib[pos])/L;
                b3 = (Ib[pos] - value[i])/L;
                break;
            }
        }

        O[3*pos  ] = checkPixelRange( a1*a2*a3*Jr[num*num*idx_r+num*idx_g+idx_b][pos] + a1*a2*b3*Jr[num*num*idx_r+num*idx_g+(idx_b+1)][pos] + a1*b2*a3*Jr[num*num*idx_r+num*(idx_g+1)+idx_b][pos] + a1*b2*b3*Jr[num*num*idx_r+num*(idx_g+1)+(idx_b+1)][pos] + b1*a2*a3*Jr[num*num*(idx_r+1)+num*idx_g+idx_b][pos] + b1*a2*b3*Jr[num*num*(idx_r+1)+num*idx_g+(idx_b+1)][pos] + b1*b2*a3*Jr[num*num*(idx_r+1)+num*(idx_g+1)+idx_b][pos] + b1*b2*b3*Jr[num*num*(idx_r+1)+num*(idx_g+1)+(idx_b+1)][pos] );
        O[3*pos+1] = checkPixelRange( a1*a2*a3*Jg[num*num*idx_r+num*idx_g+idx_b][pos] + a1*a2*b3*Jg[num*num*idx_r+num*idx_g+(idx_b+1)][pos] + a1*b2*a3*Jg[num*num*idx_r+num*(idx_g+1)+idx_b][pos] + a1*b2*b3*Jg[num*num*idx_r+num*(idx_g+1)+(idx_b+1)][pos] + b1*a2*a3*Jg[num*num*(idx_r+1)+num*idx_g+idx_b][pos] + b1*a2*b3*Jg[num*num*(idx_r+1)+num*idx_g+(idx_b+1)][pos] + b1*b2*a3*Jg[num*num*(idx_r+1)+num*(idx_g+1)+idx_b][pos] + b1*b2*b3*Jg[num*num*(idx_r+1)+num*(idx_g+1)+(idx_b+1)][pos] );
        O[3*pos+2] = checkPixelRange( a1*a2*a3*Jb[num*num*idx_r+num*idx_g+idx_b][pos] + a1*a2*b3*Jb[num*num*idx_r+num*idx_g+(idx_b+1)][pos] + a1*b2*a3*Jb[num*num*idx_r+num*(idx_g+1)+idx_b][pos] + a1*b2*b3*Jb[num*num*idx_r+num*(idx_g+1)+(idx_b+1)][pos] + b1*a2*a3*Jb[num*num*(idx_r+1)+num*idx_g+idx_b][pos] + b1*a2*b3*Jb[num*num*(idx_r+1)+num*idx_g+(idx_b+1)][pos] + b1*b2*a3*Jb[num*num*(idx_r+1)+num*(idx_g+1)+idx_b][pos] + b1*b2*b3*Jb[num*num*(idx_r+1)+num*(idx_g+1)+(idx_b+1)][pos] );
    }

    // Release memory
    for(int i = 0; i < num*num*num; ++i) {
        delete [] Jr[i];
        delete [] Jg[i];
        delete [] Jb[i];
    }
    delete [] Jr;
    delete [] Jg;
    delete [] Jb;
    delete [] Jw;
    delete [] value;
    delete [] tmp;
    delete [] tmp_r;
    delete [] tmp_g;
    delete [] tmp_b;
    delete [] tmp1;
    delete [] tmp2;
    delete [] tmp3;
    delete [] Ir;
    delete [] Ig;
    delete [] Ib;
}

void BilateralFilterC::multiChPBFICMethod_modified(uchar* I, uchar* O, int width, int height, int num)
{
    int numPixel = height*width;
	uchar* Ir = new uchar[numPixel];
	uchar* Ig = new uchar[numPixel];
	uchar* Ib = new uchar[numPixel];
	for(int i = 0; i < numPixel; ++i) {
		Ir[i] = I[3*i+0];
		Ig[i] = I[3*i+1];
		Ib[i] = I[3*i+2];
	}

    // Initialize uniform grid PBFIC
    int* value = new int[num];
    for(int i = 0; i < num; ++i) {
        value[i] = (255*i)/(num-1);
    }
    int numSample = num*num*num;
    vec3<int>* center = new vec3<int>[numSample];
    for (int i = 0; i < num; ++i) {
    for (int j = 0; j < num; ++j) {
    for (int k = 0; k < num; ++k) {
        center[num*num*i+num*j+k].x = value[i];
        center[num*num*i+num*j+k].y = value[j];
        center[num*num*i+num*j+k].z = value[k];
    }
    }
    }
    // Filtered components
    double** Jr = new double*[numSample];
    double** Jg = new double*[numSample];
    double** Jb = new double*[numSample];
    double* Jw = new double[numPixel];
    // Some tmp buffers
	double* tmp  = new double[numPixel];
    double* tmp_r  = new double[numPixel];
    double* tmp_g  = new double[numPixel];
    double* tmp_b  = new double[numPixel];
    double* tmp1 = new double[numPixel];
    double* tmp2 = new double[numPixel];
    double* tmp3 = new double[numPixel];

    // Compute bilateral filtering for each color exemplar
    for(int s = 0; s < numSample; ++s) {
        Jr[s] = new double[numPixel];
        Jg[s] = new double[numPixel];
        Jb[s] = new double[numPixel];
        for(int pos = 0; pos < numPixel; ++pos) {
            tmp[pos] = gr_kernel[abs(center[s].x-Ir[pos])]*gr_kernel[abs(center[s].y-Ig[pos])]*gr_kernel[abs(center[s].z-Ib[pos])];
        }
        //gaussian1D_recursive(tmp, Jw, width, height, sigma_s, tmp1, tmp2, tmp3);
        //cout << tmp[50] << " " << tmp[100] << " " << tmp[150] << endl;
        boxfilter(tmp, Jw, width, height, r, tmp1);
        //cout << Jw[50] << " " << Jw[100] << " " << Jw[150] << endl;
        //box_filter(tmp, Jw, width, height, r);
        //cout << Jw[50] << " " << Jw[100] << " " << Jw[150] << endl;
        //cout << "------------------------" << endl;
        for(int pos = 0; pos < numPixel; ++pos) {
            tmp_r[pos] = tmp[pos]*double(Ir[pos]);
            tmp_g[pos] = tmp[pos]*double(Ig[pos]);
            tmp_b[pos] = tmp[pos]*double(Ib[pos]);
        }
        //gaussian1D_recursive(tmp_r, Jr[s], width, height, sigma_s, tmp1, tmp2, tmp3);
        //gaussian1D_recursive(tmp_g, Jg[s], width, height, sigma_s, tmp1, tmp2, tmp3);
        //gaussian1D_recursive(tmp_b, Jb[s], width, height, sigma_s, tmp1, tmp2, tmp3);

        boxfilter(tmp_r, Jr[s], width, height, r, tmp1);
        boxfilter(tmp_g, Jg[s], width, height, r, tmp1);
        boxfilter(tmp_b, Jb[s], width, height, r, tmp1);

        for(int pos = 0; pos < numPixel; ++pos) {
            Jr[s][pos] /= Jw[pos];
            Jg[s][pos] /= Jw[pos];
            Jb[s][pos] /= Jw[pos];
        }
    }

    // Interpolation
    double* distList = new double[numSample];
    double* x_distList = new double[numSample];
    double* y_distList = new double[numSample];
    double* z_distList = new double[numSample];
    int* idList = new int[numSample];
    bool isDone;
    double accR, accG, accB, wgt, www;

    for(int pos = 0; pos < numPixel; ++pos) {
        isDone = 0;
        for (int s = 0; s < numSample; ++s) {
            if (Ir[pos] == center[s].x && Ig[pos] == center[s].y && Ib[pos] == center[s].z) {
                O[3*pos  ] = checkPixelRange(Jr[s][pos]);
                O[3*pos+1] = checkPixelRange(Jg[s][pos]);
                O[3*pos+2] = checkPixelRange(Jb[s][pos]);
                isDone = 1;
                break;
            }
        }
        if (isDone) continue;
        // compute distance
        for(int s = 0; s < numSample; ++s) {
            x_distList[s] = (double(Ir[pos])-double(center[s].x))*(double(Ir[pos])-double(center[s].x));
            y_distList[s] = (double(Ig[pos])-double(center[s].y))*(double(Ig[pos])-double(center[s].y));
            z_distList[s] = (double(Ib[pos])-double(center[s].z))*(double(Ib[pos])-double(center[s].z));
            distList[s] = x_distList[s] + y_distList[s] + z_distList[s];
            idList[s] = s;
        }

        // find the 4 nearest samples
        int nearest = 4;
        for(int k = 0; k < nearest; ++k) {
            for(int j = k+1; j < numSample; ++j) {
                if(distList[j] < distList[k]) {
                    double ddd = distList[k];
                    distList[k] = distList[j];
                    distList[j] = ddd;
                    int iii = idList[k];
                    idList[k] = idList[j];
                    idList[j] = iii;
                }
            }
        }


        accR = accG = accB = wgt = 0;
        double sigma_k = 15.0;
        for (int s = 0; s < nearest; ++s) {
            //www = 1/(distList[k]);
            www = exp(-0.5f*distList[s]/sigma_k/sigma_k);
            wgt += www;
            accR += www*Jr[idList[s]][pos];
            accG += www*Jg[idList[s]][pos];
            accB += www*Jb[idList[s]][pos];
        }

        O[3*pos  ] = checkPixelRange(accR/wgt);
        O[3*pos+1] = checkPixelRange(accG/wgt);
        O[3*pos+2] = checkPixelRange(accB/wgt);
    }


    // Release memory
    for(int i = 0; i < numSample; ++i) {
        delete [] Jr[i];
        delete [] Jg[i];
        delete [] Jb[i];
    }
    delete [] Jr;
    delete [] Jg;
    delete [] Jb;
    delete [] Jw;
    delete [] value;
    delete [] tmp;
    delete [] tmp_r;
    delete [] tmp_g;
    delete [] tmp_b;
    delete [] tmp1;
    delete [] tmp2;
    delete [] tmp3;
	delete [] Ir;
	delete [] Ig;
	delete [] Ib;
}


void BilateralFilterC::colorSampledPBFICMethod(uchar* I, uchar* O, int* RGB, int width, int height, int numSpace, int numSample)
{
	int numPixel = height*width;
	uchar* Ir = new uchar[numPixel];
	uchar* Ig = new uchar[numPixel];
	uchar* Ib = new uchar[numPixel];
    int* Ir_space = new int[numSpace];
    int* Ig_space = new int[numSpace];
    int* Ib_space = new int[numSpace];

	for(int i = 0; i < numPixel; ++i) {
		Ir[i] = I[3*i+0];
		Ig[i] = I[3*i+1];
		Ib[i] = I[3*i+2];
	}
    for(int i = 0; i < numSpace; ++i) {
       Ir_space[i] = RGB[3*i+0];
       Ig_space[i] = RGB[3*i+1];
       Ib_space[i] = RGB[3*i+2];
    }

	vec3<int>* center = new vec3<int>[numSample];
    // Fixed radius
    //PoissonDisk(Ir_space, Ig_space, Ib_space, center, numSpace, numSample, 300);
    // Varied radius
    AutoPoissonDisk(Ir_space, Ig_space, Ib_space, center, numSpace, numSample);


    // printf("centers:\n");
    // for(int k = 0; k < numSample; ++k) {
    //     printf("[%d]: %d %d %d\n", k+1, center[k].x, center[k].y, center[k].z);
    // }

    // // Visualize the clustering result
    // uchar* clusterMap = new uchar[numPixel*3];
    // vec3<double> d;
    // for(int i = 0; i < numPixel; ++i) {
    //     double min_distance = 999999;
    //     int min_k = 0;
    //     for (int k = 0; k < numSample; ++k) {
    //         d.x = Ir[i]-center[k].x;
    //         d.y = Ig[i]-center[k].y;
    //         d.z = Ib[i]-center[k].z;
    //         double distance = d.x*d.x + d.y*d.y + d.z*d.z;
    //         if (distance < min_distance) {
    //             clusterMap[3*i+0] = uchar(center[k].x);
    //             clusterMap[3*i+1] = uchar(center[k].y);
    //             clusterMap[3*i+2] = uchar(center[k].z);
    //             min_distance = distance;
    //             min_k = k;
    //         }
    //     }
    // }
    // writePPM(clusterMap, width, height, "clusterMap.ppm");

	// Filtered components
	double** Jr = new double*[numSample];
    double** Jg = new double*[numSample];
    double** Jb = new double*[numSample];
    double* Jw = new double[numPixel];
	// Some tmp buffers
    double* tmp  = new double[numPixel];
    double* tmp_r  = new double[numPixel];
    double* tmp_g  = new double[numPixel];
    double* tmp_b  = new double[numPixel];
    double* tmp1 = new double[numPixel];
    double* tmp2 = new double[numPixel];
    double* tmp3 = new double[numPixel];
	// Compute bilateral filtering for each color exemplar
    for(int k = 0; k < numSample; ++k) {
        Jr[k] = new double[numPixel];
        Jg[k] = new double[numPixel];
        Jb[k] = new double[numPixel];
        for(int i = 0; i < numPixel; ++i) {
			tmp[i] = gr_kernel[abs(center[k].x-Ir[i])]*gr_kernel[abs(center[k].y-Ig[i])]*gr_kernel[abs(center[k].z-Ib[i])];

        }
        // gaussian1D_recursive(tmp, Jw, width, height, sigma_s, tmp1, tmp2, tmp3);
        boxfilter(tmp, Jw, width, height, r, tmp1);
        for(int i = 0; i < numPixel; ++i) {
            tmp_r[i] = tmp[i]*double(Ir[i]);
            tmp_g[i] = tmp[i]*double(Ig[i]);
            tmp_b[i] = tmp[i]*double(Ib[i]);
        }
        // gaussian1D_recursive(tmp_r, Jr[k], width, height, sigma_s, tmp1, tmp2, tmp3);
        // gaussian1D_recursive(tmp_g, Jg[k], width, height, sigma_s, tmp1, tmp2, tmp3);
        // gaussian1D_recursive(tmp_b, Jb[k], width, height, sigma_s, tmp1, tmp2, tmp3);
        boxfilter(tmp_r, Jr[k], width, height, r, tmp1);
        boxfilter(tmp_g, Jg[k], width, height, r, tmp1);
        boxfilter(tmp_b, Jb[k], width, height, r, tmp1);
		for(int i = 0; i < numPixel; ++i) {
			Jr[k][i] /= Jw[i];
			Jg[k][i] /= Jw[i];
			Jb[k][i] /= Jw[i];
		}
    }

	// Interpolation
	double* distList = new double[numSample];
	double* x_distList = new double[numSample];
    double* y_distList = new double[numSample];
    double* z_distList = new double[numSample];
	int* idList = new int[numSample];
	bool isDone;
	double accR, accG, accB, wgt, www;

    for(int i = 0; i < numPixel; ++i) {
        isDone = 0;
        for (int k = 0; k < numSample; ++k) {
			if (Ir[i] == center[k].x && Ig[i] == center[k].y && Ib[i] == center[k].z) {
                O[3*i  ] = checkPixelRange(Jr[k][i]);
                O[3*i+1] = checkPixelRange(Jg[k][i]);
                O[3*i+2] = checkPixelRange(Jb[k][i]);
                isDone = 1;
				break;
            }
        }
        if (isDone) continue;
        // compute distance
        for(int k = 0; k < numSample; ++k) {
            x_distList[k] = (double(Ir[i])-double(center[k].x))*(double(Ir[i])-double(center[k].x));
			y_distList[k] = (double(Ig[i])-double(center[k].y))*(double(Ig[i])-double(center[k].y));
            z_distList[k] = (double(Ib[i])-double(center[k].z))*(double(Ib[i])-double(center[k].z));
			distList[k] = x_distList[k] + y_distList[k] + z_distList[k];
			idList[k] = k;
		}

		// find the 4 nearest samples
        int nearest = 4;
		for(int k = 0; k < nearest; ++k) {
			for(int j = k+1; j < numSample; ++j) {
				if(distList[j] < distList[k]) {
					double ddd = distList[k];
					distList[k] = distList[j];
					distList[j] = ddd;
					int iii = idList[k];
					idList[k] = idList[j];
					idList[j] = iii;
				}
			}
		}


		accR = accG = accB = wgt = 0;
		double sigma_k = 15.0;
        for (int k = 0; k < nearest; ++k) {
			//www = 1/(distList[k]);
			www = exp(-0.5f*distList[k]/sigma_k/sigma_k);
			wgt += www;
			accR += www*Jr[idList[k]][i];
			accG += www*Jg[idList[k]][i];
			accB += www*Jb[idList[k]][i];
        }

        O[3*i  ] = checkPixelRange(accR/wgt);
        O[3*i+1] = checkPixelRange(accG/wgt);
        O[3*i+2] = checkPixelRange(accB/wgt);
    }

    for(int i = 0; i < numSample; ++i) {
        delete [] Jr[i];
        delete [] Jg[i];
        delete [] Jb[i];
    }
    delete [] Jr;
    delete [] Jg;
    delete [] Jb;
    delete [] Jw;
    delete [] tmp;
    delete [] tmp_r;
    delete [] tmp_g;
    delete [] tmp_b;
    delete [] tmp1;
    delete [] tmp2;
    delete [] tmp3;
	delete [] Ir;
	delete [] Ig;
	delete [] Ib;
}
