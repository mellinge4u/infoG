/*#include "tgaimage.h"
//#include "geometry.h"
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

float z_buffer[800 * 800];

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

vector<vector<double> > matriceInverse(vector<vector<double> > matrice) {
	double det = matrice[0][0] * matrice[1][1] - matrice[0][1] * matrice[1][0];
	vector<vector<double> > matriceIA(2, vector<double>(2));
	cout << det << endl;
	matriceIA[0][0] = (matrice[1][1] / det);
	matriceIA[0][1] = (-matrice[0][1]) / det;
	matriceIA[1][0] = (-matrice[1][0]) / det;
	matriceIA[1][1] = matrice[0][0] / det;
	return matriceIA;

}

vector<vector<double> > produitMatrice2(vector<vector<double> > matrice1,
		vector<vector<double> > matrice2) {
//cout << "bonjour3 "<< endl;
	vector<vector<double> > matriceP(2, vector<double>(2));
	matriceP[0][0] = (matrice1[0][0] * matrice2[0][0])
			+ (matrice1[0][1] * matrice2[1][0]);
	matriceP[1][0] = (matrice1[1][0] * matrice2[0][0])
			+ (matrice1[1][1] * matrice2[1][0]);
	matrice1 = matriceP;
	return matriceP;
}

void triangle_filled_barycenter2(int x0, int y0, int x1, int y1, int x2, int y2,
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
//matrice A
	vector<vector<double> > matriceA(2, vector<double>(2));
	matriceA[0][0] = x1 - x0;
	matriceA[0][1] = x2 - x0;
	matriceA[1][0] = y1 - y0;
	matriceA[1][1] = y2 - y0;

//double matriceIA[][];

	vector<vector<double> > matriceIA = matriceInverse(matriceA);

//matriceB

	int xmin, xmax, ymin, ymax;
	xmin = min(min(x0, x1), x2);
	ymin = min(min(y0, y1), y2);
	xmax = max(max(x0, x1), x2);
	ymax = max(max(y0, y1), y2);
	vector<vector<double> > matriceB(2, vector<double>(1));
	int xP = 0, yP = 0;
	for (xP = xmin; xP <= xmax; xP++) {
		for (yP = ymin; yP <= ymax; yP++) {
			matriceB[0][0] = x0 - xP;
			matriceB[1][0] = y0 - yP;
			vector<vector<double> > matriceP = produitMatrice2(matriceIA,
					matriceB);

			double u = matriceP[0][0];
			//cout << "u: " << u << endl;
			double v = matriceP[1][0];
			//cout << "v: " << v << endl;
			double w = 1 - u - v;

			if (u >= 0 && v >= 0 && w >= 0) {

				cout << "bonjour4 " << endl;
				image.set(xP, yP, color);

			}
		}
	}
}

void triangle_filled_barycenter3(int x0, int y0, int z0, int x1, int y1, int z1,
		int x2, int y2, int z2, int tx0, int ty0, int tz0, int tx1, int ty1,
		int tz1, int tx2, int ty2, int tz2, TGAImage &image, TGAImage &text) {

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

	//TGAColor color = TGAColor(255 * lum, 255 * lum, 255 * lum, 255);

	// pour dessiner les triangles
	int xP, yP, txP=0, tyP = 0;
	for (xP = xmin; xP <= xmax; xP++, txP++) {
		for (yP = ymin; yP <= ymax; yP++, tyP ++) {
			vector<float> vx;
			vx.push_back(x2 - x0);
			vx.push_back(x1 - x0);
			vx.push_back(x0 - xP);

			vector<float> tvx;
			tvx.push_back(tx2 - tx0);
			tvx.push_back(tx1 - tx0);
			tvx.push_back(tx0 - txP);

			vector<float> vy;
			vy.push_back(y2 - y0);
			vy.push_back(y1 - y0);
			vy.push_back(y0 - yP);

			vector<float> tvy;
			tvy.push_back(ty2 - ty0);
			tvy.push_back(ty1 - ty0);
			tvy.push_back(ty0 - tyP);

			//produit vectoriel
			vector<float> vPv;
			vPv.push_back(vx.at(1) * vy.at(2) - vx.at(2) * vy.at(1));
			vPv.push_back(vx.at(2) * vy.at(0) - vx.at(0) * vy.at(2));
			vPv.push_back(vx.at(0) * vy.at(1) - vx.at(1) * vy.at(0));

			vector<float> tvPv;
			tvPv.push_back(tvx.at(1) * tvy.at(2) - tvx.at(2) * tvy.at(1));
			tvPv.push_back(tvx.at(2) * tvy.at(0) - tvx.at(0) * tvy.at(2));
			tvPv.push_back(tvx.at(0) * tvy.at(1) - tvx.at(1) * tvy.at(0));

			vector<float> res;
			res.push_back(1. - (vPv.at(0) + vPv.at(1)) / vPv.at(2));
			res.push_back(vPv.at(1) / vPv.at(2));
			res.push_back(vPv.at(0) / vPv.at(2));

			vector<float> tres;
			tres.push_back(1. - (tvPv.at(0) + tvPv.at(1)) / tvPv.at(2));
			tres.push_back(tvPv.at(1) / tvPv.at(2));
			tres.push_back(tvPv.at(0) / tvPv.at(2));

			if (res.at(0) >= 0 && res.at(1) >= 0 && res.at(2) >= 0) {
				float z = z0 * res.at(0) + z1 * res.at(1) + z2 * res.at(2);
				if (z > z_buffer[xP + yP * 800]) {
					z_buffer[xP + yP * 800] = z;
					TGAColor color = text.get(txP, tyP);
					color.r = color.r * lum;
					color.g = color.g * lum;
					color.b = color.b * lum;
					image.set(xP, yP, color);

				}
			}
		}
	}
}

int main(int argc, char** argv) {
	for (size_t i = 800 * 800; i--; z_buffer[i] = -numeric_limits<float>::max())
		;
	ifstream fichier(argv[1], ios::in);
	TGAImage image(800, 800, TGAImage::RGB);
	TGAImage text;
	text.read_tga_file(argv[2]);
	text.flip_vertically();

	if (fichier) {
		string ligne, chaine;
		string x, y, z, string1, string2, string3, point11, point12, point13,
				trash, string4, string5, string6;
		float x1, y1, z1, x2, y2, z2, x3, y3, z3, tx1, tx2, tx3,ty1, ty2, ty3,tz1, tz2, tz3, p1, p2, p3;
		int i, point1, point2, point3;
		i = 0;
		std::vector<string> points;
		std::vector<string> points2;
		std::vector<string> points3;
		while (std::getline(fichier, ligne)) {
			std::istringstream iss(ligne.c_str());
			if (ligne.size() > 0) {
				iss >> chaine >> x >> y >> z;

				if (chaine == "v") {
					std::istringstream iss2(x.c_str());
					std::istringstream iss3(y.c_str());
					std::istringstream iss4(z.c_str());
					iss2 >> x1;
					iss3 >> y1;
					iss4 >> z1;
					std::cout << "OK" << endl;
					//line(x,y,0,0,image,red);
					image.set(x1 * 400 + 400, y1 * 400 + 400, red);
					points.push_back(ligne);
				}

				if (chaine == "vt") {
					std::istringstream iss2(x.c_str());
					std::istringstream iss3(y.c_str());
					std::istringstream iss4(z.c_str());
					iss2 >> x1;
					iss3 >> y1;
					iss4 >> z1;
					std::cout << "OK" << endl;
					points2.push_back(ligne);
				}

				if (chaine == "f") {
					std::istringstream iss2(x.c_str());
					std::istringstream iss3(y.c_str());
					std::istringstream iss4(z.c_str());
					iss2 >> point1 >> point11;
					int pos = point11.find("/");
					int pos2 = point11.find("/", pos + 1);
					point11 = point11.substr(pos + 1, pos2 - 1);
					std::cout << "point11:" << point11 << endl;
					iss3 >> point2 >> point12;
					int pos3 = point12.find("/");
					int pos4 = point12.find("/", pos3 + 1);
					point12 = point12.substr(pos3 + 1, pos4 - 1);
					std::cout << "point12:" << point12 << endl;
					iss4 >> point3 >> point13;
					int pos5 = point13.find("/");
					int pos6 = point13.find("/", pos5 + 1);
					point13 = point13.substr(pos5 + 1, pos6 - 1);
					std::cout << "point13:" << point13 << endl;
					string1 = points[point1 - 1];
					string2 = points[point2 - 1];
					string3 = points[point3 - 1];
					p1 = stoi(point11);
					p2 = stoi(point12);
					p3 = stoi(point13);
					string4 = points2[p1 - 1];
					string5 = points2[p2 - 1];
					string6 = points2[p3 - 1];
					std::istringstream iss5(string1.c_str());
					std::istringstream iss6(string2.c_str());
					std::istringstream iss7(string3.c_str());
					iss5 >> chaine >> x1 >> y1 >> z1;
					iss6 >> chaine >> x2 >> y2 >> z2;
					iss7 >> chaine >> x3 >> y3 >> z3;
					std::cout << "point1:" << point1 << endl;
					std::cout << "point2:" << point2 << endl;
					std::cout << "point3:" << point3 << endl;
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
					std::cout << "tX1:" << tx1 << endl;
										std::cout << "tX2:" << tx2 << endl;
										std::cout << "tX3:" << tx3 << endl;
					//line(x1*400+400,y1*400+400,x2*400+400,y2*400+400,image, red);
					//line(x2*400+400,y2*400+400,x3*400+400,y3*400+400,image, red);
					//line(x1*400+400,y1*400+400,x3*400+400,y3*400+400,image, red);
					//triangle_filled_linesweep(x1 * 400 + 400, y1 * 400 + 400,
					//	x2 * 400 + 400, y2 * 400 + 400, x3 * 400 + 400,
					//	y3 * 400 + 400, image, white);
					triangle_filled_barycenter3(x1 * 400 + 400, y1 * 400 + 400,
							z1 * 400 + 400, x2 * 400 + 400, y2 * 400 + 400,
							z2 * 400 + 400, x3 * 400 + 400, y3 * 400 + 400,
							z3 * 400 + 400,tx1 * 400 + 400, ty1 * 400 + 400,
							tz1 * 400 + 400,tx2 * 400 + 400, ty2 * 400 + 400,
							tz2 * 400 + 400,tx3 * 400 + 400, ty3 * 400 + 400,
							tz3 * 400 + 400, image,text);
					std::cout << "X1:" << x1 << endl;
					std::cout << "X2:" << x2 << endl;
					std::cout << "X3:" << x3 << endl;
					//points.push_back(ligne);
				}
				i++;
				std::cout << "i:" << i << endl;
			}
		}
		std::cout << "tX1f:" << points2.at(1249) << endl;
		//line(10,20,40,50,image,red);
		// TGAImage image(800, 800, TGAImage::RGB);
		//	image.set(52, 41, red);
		image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.write_tga_file("output.tga");
		//triangle_filled_linesweep(10,30,233,2324,24234,24234,image2,red);
		//image2.flip_vertically();
		//image2.write_tga_file("triangle.tga");
		fichier.close();
	}
	return 0;
}
*/
