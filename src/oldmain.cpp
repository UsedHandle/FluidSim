#include <cstdlib>
#include <array>
#include <cmath>
#include <print>

#include <SDL3/SDL_render.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <Eigen/Dense>

constexpr int kScreenWidth{ 640 };
constexpr int kScreenHeight{ 640 };

struct Cell {
	Eigen::Vector2f U;
	float p;
};

template<size_t NX, size_t NY>
struct Cells {
	std::array<std::array<Eigen::Vector2f, NX>, NY> U;
	std::array<std::array<float, NX>, NY> p;
};

template<size_t NX, size_t NY, size_t NROWS>
struct ColomnMat {
	Eigen::Matrix<float, NX*NY, NROWS> mat;
	constexpr size_t getNX() const {
		return NX;
	}
	constexpr size_t getNY() const {
		return NY;
	}
	constexpr size_t getNROWS() const {
		return NROWS;
	}
};

const Eigen::Vector2f origin(0.f, 0.f);

const float Lx = 10.f;
const float Ly = 10.f;

const size_t nx = 50;
const size_t ny = 50;
const float dx = Lx / (float)(nx-1);
const float dy = Lx / (float)(ny-1);

const size_t nit = 50;

Eigen::Vector3f infernoCM(float t) {
	//https://observablehq.com/@flimsyhat/webgl-color-maps
	using vec3 = Eigen::Vector3f;
	const vec3 c0 = vec3(0.0002189403691192265, 0.001651004631001012, -0.01948089843709184);
	const vec3 c1 = vec3(0.1065134194856116, 0.5639564367884091, 3.932712388889277);
	const vec3 c2 = vec3(11.60249308247187, -3.972853965665698, -15.9423941062914);
	const vec3 c3 = vec3(-41.70399613139459, 17.43639888205313, 44.35414519872813);
	const vec3 c4 = vec3(77.162935699427, -33.40235894210092, -81.80730925738993);
	const vec3 c5 = vec3(-71.31942824499214, 32.62606426397723, 73.20951985803202);
	const vec3 c6 = vec3(25.13112622477341, -12.24266895238567, -23.07032500287172);

	return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

const std::array<Eigen::Vector2f, 4> normals = {
	Eigen::Vector2f(0.f,  1.f),
	Eigen::Vector2f(0.f, -1.f),
	Eigen::Vector2f(1.f,  0.f),
	Eigen::Vector2f(-1.f, 0.f)
};

const std::array<float, 4> dS = {
	dx,
	dx,
	dy,
	dy
};

// <P, N, S, W, E>
template<class T>
using Mat = std::array<T, 5>;


template<class T, size_t nx, size_t ny>
auto getMatFromIndex(
	const std::array<std::array<T, nx>,ny>& array,
	size_t i, size_t j)
{
	return Mat<T>{
		array[i][j],
		array[i+1][j],
		array[i-1][j],
		array[i][j+1],
		array[i][j-1]
	};
}

template<class T, size_t nx, size_t ny>
auto getMatFromIndexPeriodic(
	const std::array<std::array<T, nx>, ny>& array,
	size_t i, size_t j)
{
	const size_t nx = array[0].size();
	const size_t ny = array.size();
	return Mat<T>{
		array[i][j],
		array[(i+1)%ny][j],
		array[(i-1)%ny][j],
		array[i][(j + 1)%nx],
		array[i][(j - 1)%nx]
	};
}



typedef Mat<Eigen::Vector2f> UMat; // <U_P, U_N, U_S, U_W, U_E>
typedef Mat<float> PMat; // <p_P, p_N, p_S, p_W, p_E>
template<typename T>
auto interpFaceValues(const Mat<T>& mat) {
	std::array<T, 4> face;
	for (int i = 0; i < 4; ++i)
		face[i] = (mat[0] + mat[i + 1]) / 2.f;
	return face;
}

float divergence(const UMat& U) {
	using Eigen::Vector2f;
	// 1/Vol(\omega) * \sum^{faces} U_f \dot n dS
	std::array<Vector2f, 4> Uface = interpFaceValues(U);
	
	float sum = 0.f;
	for (int i = 0; i < 4; ++i)
		sum += Uface[i].dot(normals[i]) * dS[i];
	return sum/(dx*dy);
}

Eigen::Vector2f gradient(const PMat& P) {
	// 1/Vol(\omega) * \sum^{faces} P_f n dS
	std::array<float, 4> faceVal = interpFaceValues(P);
	using Eigen::Vector2f;

	Vector2f grad(0.f, 0.f);
	for (int i = 0; i < 4; ++i)
		grad += faceVal[i] * normals[i] * dS[i];
	return grad / (dx * dy);
}


Eigen::Vector2f Hterm(const UMat& U) {
	// Euler's equation so only convective term (div U dyad prod U)
	// is used because there is no diffusive term (visousity*laplace U)
	// H = div (U (x) U)
	// 1/Vol(\omega) * \sum^{faces} U_f(U_f \dot n) dS
	using Eigen::Vector2f;

	
	std::array<Vector2f, 4> Uface = interpFaceValues(U);
	
	Vector2f Hterm(0.f, 0.f);
	for (int i = 0; i < 4; ++i)
		Hterm += Uface[i] * Uface[i].dot(normals[i]) * dS[i];
	return Hterm/(dx*dy);
}

// laplace p = RHS
// solves for p
template<size_t NX, size_t NY>
std::array<std::array<float, NX>, NY> poissonEquation1it(
	const std::array<std::array<float, NX>, NY>& p,
	const std::array<std::array<Eigen::Vector2f, NX>, NY>& RHS) {
	// p = {
	//	a b c
	//  d e f
	// }
	// column mat = {
	//	a
	//	b
	//	c
	//	...
	// }

	Eigen::Matrix<float, NX * NY, NX * NY> AMat;
	for (size_t i = 0; i < ny; ++i) {
		for (size_t j = 0; j < nx; ++j) {

		}
	}
}

template<size_t NX, size_t NY>
void applyBc(Cells<NX, NY>& array) {
	using Eigen::Vector2f;
	for (size_t i = 0; i < array.size(); ++i) {
		if (i == array.size() - 1) {
			array.U[i][0] = Vector2f(1.f, 0.f);
			array.p[i][0] = 0.f;
			continue;
		}
		if (i == 0) {
			array.U[i][0] = Vector2f(0.f, 0.f);
			array.p[i][0] = array.p[i][1];
			continue;
		}
		for (size_t j = 0; j < array[0].size(); ++j) {
			if (j == 0) {
				array.U[i][j] = Vector2f(0.f, 0.f);
				array.p[i][j] = array.p[i][j+1];
				continue;
			}
			if (j == array[0].size() - 1) {
				array.U[i][j] = Vector2f(0.f, 0.f);
				array.p[i][j] = array.p[i][j - 1];
				continue;
			}
		}
	}
}

template<size_t NX, size_t NY>
void setInitialConditions(Cells<NX, NY>& cells) {
	using Eigen::Vector2f;
	for (size_t i = 0; i < NY; ++i) {
		for (size_t j = 0; j < NX; ++j) {
			cells.U[i][j] = Vector2f(0.f, 0.f);
			cells.p[i][j] = 0.f;
		}
	}
}

#include <iostream>
int main(int argc, char* args[]) {
	SDL_Window* window{};
	//SDL_Surface* screenSurface{};
	SDL_Renderer* renderer{};

	if (!SDL_Init(SDL_INIT_VIDEO)) {
		SDL_Log("SDL could not initialize! SDL error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}
	window = SDL_CreateWindow("SDL3", kScreenWidth, kScreenHeight, 0);
	if (!window) {
		SDL_Log("Window could not be created! SDL error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}

	renderer = SDL_CreateRenderer(window, NULL);
	if (!renderer) {
		SDL_Log("Renderer could not be created! SDL error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}

	
	using Eigen::Vector2f;
	std::array<std::array<Vector2f, nx>, ny> XY;
	for (size_t i = 0; i < ny; ++i) {
		for (size_t j = 0; j < nx; ++j) {
			XY[i][j] = Vector2f((float)j * (float)dx, (float)i * (float)dy) + origin;
		}
	}
	/*
	std::array<std::array<Cell, nx>, ny> cells;
	*/
	Cells<nx, ny> cells;
	setInitialConditions(cells);

	std::array<std::array<Eigen::Vector2f, nx>, ny> sample;
	std::array<std::array<float, nx>, ny> div;
	/*
	std::array<int, 9> array;
	for (int i = 0; i < 9; ++i) array[i] = i;
	auto thingy = Eigen::Map<Eigen::Matrix3i>(array.data()).transpose();
	thingy *= 2;
	std::cout << thingy << std::endl;

	array[0] += 2;
	std::cout << thingy << std::endl;
	*/


	for (int i = 0; i < ny; ++i) {
		for (int j = 0; j < nx; ++j) {
			using vec2 = Eigen::Vector2f;
			sample[i][j] = vec2(std::cos(XY[i][j](0) + 2.f * XY[i][j](1)), std::sin(XY[i][j](0)-2.f*XY[i][j](1)));
			//sample[i][j] = vec2(XY(0), XY(1));
		}
	}

	for (int i = 0; i < ny; ++i) {
		for (int j = 0; j < nx; ++j) {
			UMat mat = getMatFromIndexPeriodic(sample, i, j);
			div[i][j] = divergence(mat);
		}
	}

	//screenSurface = SDL_GetWindowSurface(window);
	bool quit{};
	SDL_Event event;
	SDL_zero(event);
	while (quit == false) {
		while (SDL_PollEvent(&event)) {
			switch (event.type) {
				case SDL_EVENT_QUIT:
					quit = true;
					break;
			}
		}
		SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);


		const float rectWidth = (float)kScreenWidth / (float)(nx - 1);
		const float rectHeight = (float)kScreenHeight / (float)(ny - 1);
		SDL_RenderClear(renderer);
		
		for (int i = 1; i < ny; ++i) {
			for (int j = 0; j < nx; ++j) {
				const int x = j;
				const int y = (ny - i) - 1;
				using vec3 = Eigen::Vector3f;

				//const float t = vec2((float)x / (float)nx, (float)y / (float)ny).norm() * 1.f/std::sqrt(2.f);
				const float t = (cells.p[y][x]+5.f)/10.f;
				vec3 col = infernoCM(t);

				//SDL_SetRenderDrawColor(renderer, (float)x/(float)nx*255.f, (float)y/(float)ny*255.f, 0, SDL_ALPHA_OPAQUE);
				SDL_SetRenderDrawColor(renderer, col(0)*255.f, col(1)*255.f, col(2)*255.f, SDL_ALPHA_OPAQUE);
				SDL_FRect rectangle = SDL_FRect{ rectWidth*(float)j, rectHeight*(float)(i-1), rectWidth, rectHeight};
				SDL_RenderFillRect(renderer, &rectangle);
			}
		}
		SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
		for (int i = 0; i < ny; ++i) {
			for (int j = 0; j < nx; ++j) {
				const int x = j;
				const int y = (ny - i) - 1;
				using vec2 = Eigen::Vector2f;
				const Eigen::Matrix2f transpose{
					{ 1.f,  0.f },
					{ 0.f, -1.f }
				};
				{
					const vec2 F = transpose * (cells.U[y][x]);
					const vec2 origin(rectWidth * (float)j, rectHeight * (float)(i));
					const vec2 end = origin + F * 5.f;
					SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
					SDL_RenderLine(renderer, origin(0), origin(1), end(0), end(1));
					SDL_SetRenderDrawColor(renderer, 0, 255, 0, SDL_ALPHA_OPAQUE);
					SDL_RenderPoint(renderer, origin(0), origin(1));
				}
				
				if ((i + j) % 2== 0) {
					UMat mat = getMatFromIndexPeriodic(cells.U, y, x);
					std::array<vec2, 4> faces = interpFaceValues(mat);
					for(int facei = 0; facei < 4; ++facei){
						vec2 F = transpose * faces[facei];
						const vec2 faceOrigin(rectWidth* ((float)j + normals[facei](0)*0.5f), rectHeight* ((float)i + normals[facei](1)*0.5f));
						const vec2 end = faceOrigin + F * 5.f;
						SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE);
						SDL_RenderLine(renderer, faceOrigin(0), faceOrigin(1), end(0), end(1));
						SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE);
						SDL_RenderPoint(renderer, faceOrigin(0), faceOrigin(1));
					}

				}
			}
		}
		SDL_RenderPresent(renderer);
	}
	SDL_DestroyRenderer(renderer);
	//SDL_DestroySurface(screenSurface);
	SDL_DestroyWindow(window);
	SDL_Quit();
}