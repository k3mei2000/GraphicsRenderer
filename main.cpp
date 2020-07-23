#include <vector>
#include <cmath>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
const int width = 800;
const int height = 800;

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

void line (Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
	bool isSteep = false;
	if (std::abs(p0.x-p1.x) < std::abs(p0.y-p1.y)) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		isSteep = true;
	}

	if (p0.x > p1.x) {
		std::swap(p0,p1);
	}

	int dx = p1.x - p0.x;
	int dy = p1.y - p0.y;

	int derror = std::abs(dy) * 2;
	int error = 0;
	
	int y = p0.y;
	for (int x = p0.x; x<=p1.x; x++) {
		if (isSteep) {
			image.set(y,x,color);
		} else {
			image.set(x,y,color);
		}

		error += derror;
		if (error > dx) {
			y += (p1.y>p0.y?1:-1);
			error -= dx * 2;
		}
	}
}
// Triangle Method Using Line Sweeping
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
	if (t0.y > t1.y) std::swap(t0,t1);
	if (t0.y > t2.y) std::swap(t0,t2);
	if (t1.y > t2.y) std::swap(t1,t2);

	line(t0, t1, image, color);
	line(t1, t2, image, color);
	line(t2, t0, image, color);

	
	for (int i = t0.y; i <= t2.y; i++) {
		double alpha, beta;
		if (t2.y - t0.y == 0) {
			alpha = 0;
		} else {
			alpha = (i - t0.y) / (float) (t2.y - t0.y);
		}
	
		if (i <= t1.y) {
			if (t1.y - t0.y == 0) {
				beta = 0;
			} else {
				beta = (i - t0.y) / (float) (t1.y - t0.y);
			}
		} else {
			if (t2.y - t1.y == 0) {
				beta = (i - t1.y) / (float) (t2.y - t1.y);
			} else {
				beta = (i - t1.y) / (float) (t2.y - t1.y);
			}
			
		}
		
		int leftX, rightX;
		if (t0.x <= t2.x) {
			leftX = t0.x + alpha * (t2.x - t0.x);
		} else {
			leftX = t0.x - alpha * (t0.x - t2.x);
		}

		if (i <= t1.y) {
			if (t0.x <= t1.x) {
				rightX = t0.x + beta * (t1.x - t0.x);
			} else {
				rightX = t0.x - beta * (t0.x - t1.x);
			}
		} else {
			if (t1.x <= t2.x) {
				rightX = t1.x + beta * (t2.x - t1.x);
			} else {
				rightX = t1.x - beta * (t1.x - t2.x);
			}
		}

		if (leftX > rightX) {
			std::swap(leftX, rightX);
		}

		for (int j = leftX; j <= rightX; j++) {
			Vec2i leftPoint = Vec2i(leftX, i);
			Vec2i rightPoint = Vec2i(rightX, i);
			line(leftPoint,rightPoint,image,color);
		}
	}
}
/*
Vec3f barycentric(Vec2i *pts, Vec2i P) {
	Vec3f u = cross(Vec3f(pts[1][0] - pts[0][0], pts[2][0] - pts[0][0], pts[0][0] - P[0]),
		Vec3f(pts[1][1] - pts[0][1], pts[2][1] - pts[0][1], pts[0][1] - P[1]));
	
	if (abs(u[2]) < 1) return Vec3f(-1,1,1);
	return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
}

// Triangle method using barycentric coordinates.
void triangle(Vec2i *pts, TGAImage &image, TGAColor color) {
	Vec2i bboxmin(image.get_width()-1, image.get_height()-1);
	Vec2i bboxmax(0,0);
	Vec2i clamp(image.get_width()-1, image.get_height() - 1);
	for (int i=0; i<3; i++) {
		for (int j=0; j<2; j++) {
			bboxmin[j] = std::max(0, std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
		}
	}
	Vec2i P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts, P);
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) {
				image.set(P.x,P.y, color);
			}
		}
	}
}
*/

Vec3f barycentric(Vec3f *pts, Vec3f P) {
	Vec3f s[2];
	for (int i = 2; i--; ) {
		s[i][0] = pts[1][i] - pts[0][i];
		s[i][1] = pts[2][i] - pts[0][i];
		s[i][2] = pts[0][i] - P[i];
	}
	Vec3f u = cross(s[0],s[1]);
	if (abs(u[2]) > 1e-2) return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
	return Vec3f(-1,1,1);
}

void triangle(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color) {
	Vec2f bboxmin( std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
	Vec2f clamp(image.get_width()-1, image.get_height()-1);

	for (int i=0; i<3; i++) {
		for (int j=0; j<2; j++) {
			bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
		}
	}
	Vec3f P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts, P);
			if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
			P.z =0;
			for (int i = 0; i<3; i++) P.z += pts[i][2]*bc_screen[i];
			if (zbuffer[int(P.x+P.y*width)]<P.z) {
				zbuffer[int(P.x+P.y*width)] = P.z;
				image.set(P.x,P.y,color);
			}
		}
	}
}

void triangle(Vec3f *pts, Vec3f *texture_pts, float intensity, float *zbuffer, TGAImage &image, TGAImage &textureImage) {
	Vec2f bboxmin( std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
	Vec2f clamp(image.get_width()-1, image.get_height()-1);
	float texture_width = textureImage.get_width();
	float texture_height = textureImage.get_height();

	for (int i=0; i<3; i++) {
		for (int j=0; j<2; j++) {
			bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
		}
	}
	Vec3f P;
	Vec3f texture;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts, P);
			if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
			P.z =0;
			texture.x = 0;
			texture.y = 0;
			for (int i = 0; i<3; i++) {
				P.z += pts[i][2]*bc_screen[i];
				texture.x += texture_pts[(i)%3][0]*bc_screen[(3-i)%3];
				texture.y += texture_pts[(i)%3][1]*bc_screen[(3-i)%3];
			}
			

			if (zbuffer[int(P.x+P.y*width)]<P.z) {
				zbuffer[int(P.x+P.y*width)] = P.z;
				TGAColor color = textureImage.get((int) (texture.x * texture_width),(int) (texture.y * texture_height));
				image.set(P.x,P.y,color);
			}
		}
	}
}


Vec3f world2screen(Vec3f v) {
	return Vec3f(int((v.x+1.)*width/2. + .5), int((v.y+1.)*height/2.+.5), v.z);
}

int main(int argc, char** argv) {
	
	if (2==argc) {
		model = new Model(argv[1]);
	} else {
		model = new Model("obj/african_head.obj");
	}
	
	bool vt_success = false;
	TGAImage texture = TGAImage();
	if (texture.read_tga_file("obj/african_head_diffuse.tga")) {
		vt_success = true;
		texture.flip_vertically();
	}

	TGAImage image(width, height, TGAImage::RGB);
	Vec3f light_dir(0,0,-1);
	float *zbuffer = new float[width*height];
	for (int i =width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
	/*
	line(0,0,20,50,image, white);
	line(0,0,50,20,image, red);
	image.set(52, 41, red);
	*/
	/*
	for (int i=0; i<model->nfaces(); i++) {
		std::vector<int> face = model->face(i);
		for (int j=0; j<3; j++) {
			Vec3f v0 = model->vert(face[j]);
			Vec3f v1 = model->vert(face[(j+1)%3]);
			int x0 = (v0.x+1.)*width/2.;
			int y0 = (v0.y+1.)*height/2.;
			int x1 = (v1.x+1.)*width/2.;
			int y1 = (v1.y+1.)*height/2.;
			line(x0, y0, x1, y1, image, white);
		}
	}
	*/
	
	
	for (int i = 0; i < model->nfaces(); i++) {
		
		std::vector<int> face = model->face(i);
		Vec3f pts[3], world_coords[3], texture_pts[3];
		for (int j=0; j<3; j++) {
			world_coords[j] = model->vert(face[2*j]);
			pts[j] = world2screen(world_coords[j]);

			texture_pts[j] = model->texture(face[2*j+1]);
		}


		
		Vec3f n = cross((world_coords[2]-world_coords[0]),(world_coords[1]-world_coords[0]));
		n.normalize();
		float intensity = n * light_dir;
		if (intensity > 0) {
			if (vt_success) {
				triangle(pts, texture_pts, intensity, zbuffer, image, texture);
			} else {
				float brightness = intensity * 255;
				triangle(pts, zbuffer, image, TGAColor(brightness,brightness,brightness,255));
			}
			
		}
		
	}
	

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	delete model;
	return 0;
}

