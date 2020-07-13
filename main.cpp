#include "tgaimage.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
	bool isSteep = false;
	if (std::abs(x0-x1) < std::abs(y0-y1)) {
		std::swap(x0,y0);
		std::swap(x1,y1);
		isSteep = true;
	}

	if (x0 > x1) {
		std::swap(x0,x1);
		std::swap(y0,y1);
	}

	int dx = x1 - x0;
	int dy = y1 - y0;

	int derror = std::abs(dy) * 2;
	int error = 0;
	
	int y = y0;
	for (int x = x0; x<=x1; x++) {
		if (isSteep) {
			image.set(y,x,color);
		} else {
			image.set(x,y,color);
		}

		error += derror;
		if (error > dx) {
			y += (y1>y0?1:-1);
			error -= dx * 2;
		}
	}
}

int main(int argc, char** argv) {
	TGAImage image(100, 100, TGAImage::RGB);
	line(0,0,20,50,image, white);
	line(0,0,50,20,image, red);
	image.set(52, 41, red);

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}

