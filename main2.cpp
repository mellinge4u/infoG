#include "tgaimage.h"
#include <fstream>
#include <string>
#include <fstream>
#include <math.h>
#include <complex>
#include <iostream>
#include <istream>
#include <vector>
#include <sstream>
#include <limits>

using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
const int HEIGHT = 800;
const int WIDTH = 800;

vector<string> f_points;
vector<string> f_points2;

float cam[3] = { 0, 0, 10 };

float z_buffer[800 * 800];

void f_points_display() {
	for (int i = 0; i < f_points.size(); i++) {
		cout << "f_points: " << i << f_points[i] << endl;
	}
}

void f_points2_display() {
	for (int i = 0; i < f_points.size(); i++) {
		cout << "f_points2: " << i << f_points2[i] << endl;
	}
}

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
	bool steep = false;
	if (abs(x0 - x1) < abs(y0 - y1)) {
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if (x0 > x1) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}
	int dx = x1 - x0;
	int dy = y1 - y0;
	int derror2 = abs(dy) * 2;
	int error2 = 0;
	int y = y0;
	for (int x = x0; x <= x1; x++) {
		if (steep) {
			image.set(y, x, color);
		} else {
			image.set(x, y, color);
		}
		error2 += derror2;
		if (error2 > dx) {
			y += (y1 > y0 ? 1 : -1);
			error2 -= dx * 2;
		}
	}
}

void triangle_filled_linesweep(int x0, int y0, int x1, int y1, int x2, int y2,
		TGAImage &image, TGAColor color) {

	if (y0 > y1) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}
	if (y0 > y2) {
		std::swap(x0, x2);
		std::swap(y0, y2);
	}
	if (y1 > y2) {
		std::swap(x1, x2);
		std::swap(y1, y2);
	}
	int hauteur_max = y2 - y0;
	if (hauteur_max != 0) {
		for (int i = y0; i <= y1; i++) {
			int hauteur_segment = y1 - y0 + 1;
			if (hauteur_segment != 0) {
				float a = (float) (i - y0) / hauteur_max;
				float b = (float) (i - y0) / hauteur_segment;
				int xA = x0 + (x2 - x0) * a;
				int yA = y0 + (y2 - y0) * a;
				int xB = x0 + (x1 - x0) * b;
				int yB = y0 + (y1 - y0) * b;
				if (xA > xB) {
					std::swap(xA, xB);
					std::swap(yA, yB);
				}
				for (int j = xA; j <= xB; j++) {
					image.set(j, i, color);
				}
			}
		}
		for (int i = y1; i <= y2; i++) {
			int hauteur_segment = y2 - y1 + 1;
			if (hauteur_segment != 0) {
				float alpha = (float) (i - y0) / hauteur_max;
				float beta = (float) (i - y1) / hauteur_segment;
				int xA = x0 + (x2 - x0) * alpha;
				int yA = y0 + (y2 - y0) * alpha;
				int xB = x1 + (x2 - x1) * beta;
				int yB = y1 + (y2 - y1) * beta;
				if (xA > xB) {
					std::swap(xA, xB);
					std::swap(yA, yB);
				}
				for (int j = xA; j <= xB; j++) {
					image.set(j, i, color);
				}
			}
		}
	}
}

void triangle_filled_bary(vector<int> v1, vector<int> v2, vector<int> text_v,
		TGAImage &image, TGAColor color) {
	int x0 = v1[0];
	int y0 = v1[1];
	int z0 = v1[2];

	int x1 = v2[0];
	int y1 = v2[1];
	int z1 = v2[2];

	int x2 = text_v[0];
	int y2 = text_v[1];
	int z2 = text_v[2];

	int xmin, xmax, ymin, ymax;
	xmin = min(min(x0, x1), x2);
	ymin = min(min(y0, y1), y2);
	xmax = max(max(x0, x1), x2);
	ymax = max(max(y0, y1), y2);
	// pour l'éclairage
	vector<float> vx;
	vx.push_back(x1 - x0);
	vx.push_back(y1 - y0);
	vx.push_back(z1 - z0);

	vector<float> vy;
	vy.push_back(x2 - x0);
	vy.push_back(y2 - y0);
	vy.push_back(z2 - z0);

	vector<float> vPv;
	vPv.push_back(vx.at(1) * vy.at(2) - vx.at(2) * vy.at(1));
	vPv.push_back(vx.at(2) * vy.at(0) - vx.at(0) * vy.at(2));
	vPv.push_back(vx.at(0) * vy.at(1) - vx.at(1) * vy.at(0));

	//normalisation
	float normalisation = sqrtf(
			pow(vPv.at(0), 2) + pow(vPv.at(1), 2) + pow(vPv.at(2), 2));
	vPv.at(0) /= normalisation;
	vPv.at(1) /= normalisation;
	vPv.at(2) /= normalisation;

	//produit scalaire
	float lum = abs(vPv.at(0) * 0 + vPv.at(1) * 0 + vPv.at(2) * 1);

	color = TGAColor(255 * lum, 255 * lum, 255 * lum, 255);

	// pour dessiner les triangles
	int xP, yP, txP = 0, tyP = 0;
	for (xP = xmin; xP <= xmax; xP++, txP++) {
		for (yP = ymin; yP <= ymax; yP++, tyP++) {
			vector<float> vx;
			vx.push_back(x2 - x0);
			vx.push_back(x1 - x0);
			vx.push_back(x0 - xP);

			vector<float> vy;
			vy.push_back(y2 - y0);
			vy.push_back(y1 - y0);
			vy.push_back(y0 - yP);

			vector<float> tvy;

			//produit vectoriel
			vector<float> vPv;
			vPv.push_back(vx.at(1) * vy.at(2) - vx.at(2) * vy.at(1));
			vPv.push_back(vx.at(2) * vy.at(0) - vx.at(0) * vy.at(2));
			vPv.push_back(vx.at(0) * vy.at(1) - vx.at(1) * vy.at(0));

			vector<float> res;
			res.push_back(1. - (vPv.at(0) + vPv.at(1)) / vPv.at(2));
			res.push_back(vPv.at(1) / vPv.at(2));
			res.push_back(vPv.at(0) / vPv.at(2));

			if (res.at(0) >= 0 && res.at(1) >= 0 && res.at(2) >= 0) {
				float z = z0 * res.at(0) + z1 * res.at(1) + z2 * res.at(2);
				if (z > z_buffer[xP + yP * 800]) {
					z_buffer[xP + yP * 800] = z;
					image.set(xP, yP, color);

				}
			}
		}
	}
}

void triangle_filled_bary2(vector<float> v1, vector<float> v2, vector<float> v3,
		vector<float> text_v1, vector<float> text_v2, vector<float> text_v3,
		TGAImage &image, TGAImage &text) {
	float x0 = v1[0];
	float y0 = v1[1];
	float z0 = v1[2];

	float x1 = v2[0];
	float y1 = v2[1];
	float z1 = v2[2];

	float x2 = v3[0];
	float y2 = v3[1];
	float z2 = v3[2];

	float tx0 = text_v1[0];
	float ty0 = text_v1[1];
	float tz0 = text_v1[2];

	float tx1 = text_v2[0];
	float ty1 = text_v2[1];
	float tz1 = text_v2[2];

	float tx2 = text_v3[0];
	float ty2 = text_v3[1];
	float tz2 = text_v3[2];

	float xmin, xmax, ymin, ymax;
	xmin = min(min(x0, x1), x2);
	ymin = min(min(y0, y1), y2);
	xmax = max(max(x0, x1), x2);
	ymax = max(max(y0, y1), y2);
	// pour l'éclairage
	vector<float> vx;
	vx.push_back(x1 - x0);
	vx.push_back(y1 - y0);
	vx.push_back(z1 - z0);

	vector<float> vy;
	vy.push_back(x2 - x0);
	vy.push_back(y2 - y0);
	vy.push_back(z2 - z0);

	vector<float> vPv;
	vPv.push_back(vx.at(1) * vy.at(2) - vx.at(2) * vy.at(1));
	vPv.push_back(vx.at(2) * vy.at(0) - vx.at(0) * vy.at(2));
	vPv.push_back(vx.at(0) * vy.at(1) - vx.at(1) * vy.at(0));

	//normalisation
	float normalisation = sqrtf(
			pow(vPv.at(0), 2) + pow(vPv.at(1), 2) + pow(vPv.at(2), 2));
	vPv.at(0) /= normalisation;
	vPv.at(1) /= normalisation;
	vPv.at(2) /= normalisation;

	//produit scalaire
	float lum = abs(vPv.at(0) * 0 + vPv.at(1) * 0 + vPv.at(2) * 1);

	//color = TGAColor(255 * lum, 255 * lum, 255 * lum, 255);

	// pour dessiner les triangles
	int xP, yP, txP = 0, tyP = 0;
	for (xP = xmin; xP <= xmax; xP++) {
		for (yP = ymin; yP <= ymax; yP++) {
			vector<float> vx;
			vx.push_back(x2 - x0);
			vx.push_back(x1 - x0);
			vx.push_back(x0 - xP);

			vector<float> vy;
			vy.push_back(y2 - y0);
			vy.push_back(y1 - y0);
			vy.push_back(y0 - yP);

			//produit vectoriel
			vector<float> vPv;
			vPv.push_back(vx.at(1) * vy.at(2) - vx.at(2) * vy.at(1));
			vPv.push_back((vx.at(2) * vy.at(0) - vx.at(0) * vy.at(2)));
			vPv.push_back(vx.at(0) * vy.at(1) - vx.at(1) * vy.at(0));

			vector<float> res;
			res.push_back(1. - (vPv.at(0) + vPv.at(1)) / vPv.at(2));
			res.push_back(vPv.at(1) / vPv.at(2));
			res.push_back(vPv.at(0) / vPv.at(2));

			// on applique le "poid" sur les coordonnées des textures
			float p_tv1x = tx0 * res[0];
			float p_tv1y = ty0 * res[0];
			float p_tv2x = tx1 * res[1];
			float p_tv2y = ty1 * res[1];
			float p_tv3x = tx2 * res[2];
			float p_tv3y = ty2 * res[2];

			if (res.at(0) >= 0 && res.at(1) >= 0 && res.at(2) >= 0) {
				float z = z0 * res.at(0) + z1 * res.at(1) + z2 * res.at(2);
				if (z > z_buffer[xP + yP * 800]) {
					z_buffer[xP + yP * 800] = z;
					float t_x = p_tv1x + p_tv2x + p_tv3x;
					float t_y = p_tv1y + p_tv2y + p_tv3y;
					TGAColor color = text.get(t_x, t_y);
					color.r = color.r * lum;
					color.g = color.g * lum;
					color.b = color.b * lum;
					image.set(xP, yP, color);
				}
			}
		}
	}
}

void mult_mat_44_44(float* m1, float* m2, float* m_res) {
	for (int i = 0; i < 16; i++)
		m_res[i] = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				m_res[i * 4 + j] += m2[j + k * 4] * m1[i * 4 + k];
			}
		}
	}
}

void mult_mat_44_41(float* m_tra, float* m_coord, float* m_res) {
	for (int i = 0; i < 4; i++)
		m_res[i] = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m_res[i] += m_tra[i * 4 + j] * m_coord[j];
		}
	}
}

void calculateTransformation(float* m_coord, float* m_res) {
	m_res[0] = m_coord[0] / m_coord[3];
	m_res[1] = m_coord[1] / m_coord[3];
	m_res[2] = m_coord[2] / m_coord[3];
}

void set_viewport(int x, int y, int w, int h, float* viewport) {
	viewport[0 * 4 + 3] = w / 2.f;
	viewport[1 * 4 + 3] = h / 2.f;
	viewport[2 * 4 + 3] = w / 2.f;
	viewport[0 * 4 + 0] = w / 2.f;
	viewport[1 * 4 + 1] = h / 2.f;
	viewport[2 * 4 + 2] = w / 2.f;
}

void applyTransformation(vector<float> &v1, vector<float> &v2,
		vector<float> &v3) {
	float m_v1[4] = { v1[0], v1[1], v1[2], 1.f };
	float m_v2[4] = { v2[0], v2[1], v2[2], 1.f };
	float m_v3[4] = { v3[0], v3[1], v3[2], 1.f };

	float projection[16] = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f,
			1.f, 0.f, 0.f, 0.f, 0.f, 1.f };
	projection[3 * 4 + 2] = -1.f / 90.;

	float viewport[16] = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f,
			1.f, 0.f, 0.f, 0.f, 0.f, 1.f };

	// a faire
	float rotation[16] = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f,
			1.f, 0.f, 0.f, 0.f, 0.f, 1.f };

	float transfo[16] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };

	set_viewport(WIDTH, HEIGHT, WIDTH, HEIGHT, viewport);

	mult_mat_44_44(viewport, projection, transfo);

	float m_inter[4] = { 0.f, 0.f, 0.f, 0.f };
	float m_f[3] = { 0.f, 0.f, 0.f };

	mult_mat_44_41(transfo, m_v1, m_inter);

	calculateTransformation(m_inter, m_f);
	v1[0] = m_f[0];
	v1[1] = m_f[1];
	v1[2] = m_f[2];

	mult_mat_44_41(transfo, m_v2, m_inter);

	calculateTransformation(m_inter, m_f);
	v2[0] = m_f[0];
	v2[1] = m_f[1];
	v2[2] = m_f[2];

	mult_mat_44_41(transfo, m_v3, m_inter);

	calculateTransformation(m_inter, m_f);
	v3[0] = m_f[0];
	v3[1] = m_f[1];
	v3[2] = m_f[2];

}

int main(int argc, char** argv) {
	for (size_t i = 800 * 800; i--; z_buffer[i] = -numeric_limits<float>::max())
		;
	ifstream fichier(argv[1], ios::in);
	TGAImage image(800, 800, TGAImage::RGB);
	TGAImage text;
	text.read_tga_file("african_head_diffuse.tga");
	int h = text.get_height();
	int w = text.get_width();
	text.flip_vertically();

	if (fichier) {
		string ligne, chaine;
		string x, y, z, string1, string2, string3, point11, point12, point13,
				trash, string4, string5, string6;
		float x1, y1, z1, x2, y2, z2, x3, y3, z3, tx1, tx2, tx3, ty1, ty2, ty3,
				tz1, tz2, tz3, p1, p2, p3;
		int i, point1, point2, point3;
		i = 0;
		while (std::getline(fichier, ligne)) {
			std::istringstream iss(ligne.c_str());
			if (ligne.size() > 0) {
				iss >> chaine >> x >> y >> z;

				if (chaine == "v") {
					f_points.push_back(ligne);
				}

				if (chaine == "vt") {
					f_points2.push_back(ligne);
				}

				if (chaine == "f") {
					std::istringstream iss2(x.c_str());
					std::istringstream iss3(y.c_str());
					std::istringstream iss4(z.c_str());
					iss2 >> point1 >> point11;
					int pos = point11.find("/");
					int pos2 = point11.find("/", pos + 1);
					point11 = point11.substr(pos + 1, pos2 - 1);
					iss3 >> point2 >> point12;
					int pos3 = point12.find("/");
					int pos4 = point12.find("/", pos3 + 1);
					point12 = point12.substr(pos3 + 1, pos4 - 1);
					iss4 >> point3 >> point13;
					int pos5 = point13.find("/");
					int pos6 = point13.find("/", pos5 + 1);
					point13 = point13.substr(pos5 + 1, pos6 - 1);
					string1 = f_points[point1 - 1];
					string2 = f_points[point2 - 1];
					string3 = f_points[point3 - 1];
					p1 = stoi(point11);
					p2 = stoi(point12);
					p3 = stoi(point13);
					string4 = f_points2[p1 - 1];
					string5 = f_points2[p2 - 1];
					string6 = f_points2[p3 - 1];
					std::istringstream iss5(string1.c_str());
					std::istringstream iss6(string2.c_str());
					std::istringstream iss7(string3.c_str());
					iss5 >> chaine >> x1 >> y1 >> z1;
					iss6 >> chaine >> x2 >> y2 >> z2;
					iss7 >> chaine >> x3 >> y3 >> z3;
					std::istringstream iss8(string1.c_str());
					std::istringstream iss9(string2.c_str());
					std::istringstream iss10(string3.c_str());
					iss8 >> chaine >> x1 >> y1 >> z1;
					iss9 >> chaine >> x2 >> y2 >> z2;
					iss10 >> chaine >> x3 >> y3 >> z3;
					std::istringstream iss11(string4.c_str());
					std::istringstream iss12(string5.c_str());
					std::istringstream iss13(string6.c_str());
					iss11 >> chaine >> tx1 >> ty1 >> tz1;
					iss12 >> chaine >> tx2 >> ty2 >> tz2;
					iss13 >> chaine >> tx3 >> ty3 >> tz3;
					vector<float> v1, v2, v3, v4, v5, v6;
					v1.push_back(x1);
					v1.push_back(y1);
					v1.push_back(z1);
					v2.push_back(x2);
					v2.push_back(y2);
					v2.push_back(z2);
					v3.push_back(x3);
					v3.push_back(y3);
					v3.push_back(z3);

					v4.push_back(tx1 * w);
					v4.push_back(ty1 * h);
					//v4.push_back(tz1 * 512 + 512);

					v5.push_back(tx2 * w);
					v5.push_back(ty2 * h);
					//v5.push_back(tz2 * 512 + 512);

					v6.push_back(tx3 * w);
					v6.push_back(ty3 * h);
					//v6.push_back(tz3 * 512 + 512);
					//triangle_filled_bary(v1, v2, v3, image, white);
					cout << " v1 x 1 " << v1[0] << endl;
					cout << " v2 x 1 " << v1[1] << endl;
					cout << " v2 x 1 " << v1[2] << endl;
					applyTransformation(v1, v2, v3);
					cout << " v1 x 2 " << v1[0] << endl;
					cout << " v2 x 2 " << v1[1] << endl;
					cout << " v2 x 2 " << v1[2] << endl;
					triangle_filled_bary2(v1, v2, v3, v4, v5, v6, image, text);
				}
				i++;
			}
		}
		//f_points2_display();
		image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.write_tga_file("output.tga");
		fichier.close();
	}
	return 0;
}
