#include "fileIO.h"

uchar checkPixelRange(double x) {
	if(x > 255)
		return uchar(255);
	else if(x < 0)
		return uchar(0);
	else
		return uchar(x);
}

void saveFloat(float* data, const int width, const int height, string filename) {

	ofstream ofs(filename.c_str(), ios::out);
	ofs << fixed;
	for(int y = 0; y < height; ++y) {
		for(int x = 0; x < width; ++x) {
			ofs << setprecision(5) << data[y*width+x] << " ";
		}
		ofs << endl;
	}
	ofs.close();

}

// Get size of pbm/pgm/ppm format
void sizePGM(int& width, int& height, string filename) {

	char ch;
	ifstream ifs(filename.c_str(), ios::binary);
	ch = ifs.get();		// eat "P"
	ch = ifs.get();		// eat "5"
	ifs >> width
		>> height;
	ifs.close();

}

void readPGM(uchar* I, string filename) {

	char ch;
	int width, height, tmp;
	ifstream ifs(filename.c_str(), ios::binary);
	if(!ifs) { cerr << "No such file..." << endl; return; }
	ch = ifs.get();		// eat "P"
	ch = ifs.get();		// eat "5"
	ifs >> width
		>> height
		>> tmp;		// tmp will be 255
	if(tmp != 255) {
		cerr << "Input range is wrong!\n";
		return;
	}
	ch = ifs.get();		// eat "\n"
	ifs.read(reinterpret_cast<char*>(I), height*width);
	ifs.close();

}

void readPPM(uchar* I, string filename) {

	char ch;
	int width, height, tmp;
	ifstream ifs(filename.c_str(), ios::binary);
	if(!ifs) { cerr << "No such file..." << endl; return; }
	ch = ifs.get();		// eat "P"
	ch = ifs.get();		// eat "6"
	ifs >> width
		>> height
		>> tmp;		// tmp will be 255
	if(tmp != 255) {
		cerr << "Input range is wrong!\n";
		return;
	}
	ch = ifs.get();		// eat "\n"
	ifs.read(reinterpret_cast<char*>(I), height*width*3);
	ifs.close();

}

void readPGM2f(float* I, string filename) {

	char ch;
	int width, height, tmp;
	ifstream ifs(filename.c_str(), ios::binary);
	if(!ifs) { cerr << "No such file..." << endl; return; }
	ch = ifs.get();		// eat "P"
	ch = ifs.get();		// eat "5"
	ifs >> width
		>> height
		>> tmp;		// tmp will be 255
	if(tmp != 255) {
		cerr << "Input range is wrong!\n";
		return;
	}
	ch = ifs.get();		// eat "\n"

	uchar* buffer = new uchar[width*height];
	ifs.read(reinterpret_cast<char*>(buffer), height*width);
	ifs.close();
	for(int i = 0; i < width*height; ++i) {
		I[i] = float(buffer[i]);
	}
	delete [] buffer;
}

void readPPM2RGBf(float* R, float* G, float* B, string filename) {

	uchar ch;
	int width, height, tmp;
	ifstream ifs(filename.c_str(), ios::binary);
	if(!ifs) { cerr << "No such file..." << endl; return; }
	ch = ifs.get();		// eat "P"
	ch = ifs.get();		// eat "6"
	ifs >> width
		>> height
		>> tmp;		// tmp will be 255
	if(tmp != 255) {
		cerr << "Input range is wrong!\n";
		return;
	}
	ch = ifs.get();		// eat "\n"

	for(int i = 0; i < width*height; ++i) {
		ch = ifs.get(); R[i] = float(ch);
		ch = ifs.get(); G[i] = float(ch);
		ch = ifs.get(); B[i] = float(ch);
	}

	ifs.close();
}

void writeRGBf2PPM(float* R, float* G, float* B, const int width, const int height, string filename)
{
	ofstream ofs(filename.c_str(), ios::binary);
	ofs << "P6" << endl
		<< width << " " << height << endl
		<< "255" << endl;
	for(int i = 0; i < width*height; ++i) {
		ofs << uchar(R[i]) << uchar(G[i]) << uchar(B[i]);
	}
	ofs.close();
}


void sizeColorSpace(int& numSpace, string filename)
{
	ifstream ifs(filename.c_str(), ios::binary);
	if(!ifs) { cerr << "No such file..." << filename << endl; return; }
	// get length of file
    int number_of_lines = 0;
    string line;
    while (getline(ifs, line))
        ++number_of_lines;
    numSpace = number_of_lines;
}

void readColorSpace(int* RGB, int numSpace, string filename)
{
	FILE *file = fopen(filename.c_str(), "r");
	for (int i = 0; i < numSpace; ++i)
		// opencv write (B,G,R)
		fscanf(file, "%d %d %d\n", &RGB[3*i+2], &RGB[3*i+1], &RGB[3*i]);
}
