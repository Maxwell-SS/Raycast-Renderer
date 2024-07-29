// SDL
#include <SDL.h>

// std
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

const int WIDTH = 640;
const int HEIGHT = 480;
const int MAPWIDTH = 24;
const int MAPHEIGHT = 24;

int m[MAPWIDTH][MAPHEIGHT] = {
	{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
	{1,0,0,0,0,0,1,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,1,0,0,0,0,0,2,0,3,0,0,0,0,0,0,3,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,1,0,0,0,0,0,2,5,2,5,2,5,2,5,2,5,0,1},
	{1,0,0,0,0,0,1,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,1},
	{1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,2,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,2,2,2,0,2,2,2,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,3,0,3,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,1},
	{1,3,3,0,3,3,3,3,3,2,2,4,3,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,3,0,3,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,0,1},
	{1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,4,0,1},
	{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,4,0,1},
	{1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,4,0,1},
	{1,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,4,4,4,4,0,1},
	{1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
};

struct vec2 {
	double x{}, y{}; // the curly brackets set the initial value to 0

	vec2() {}
	vec2(double xx, double yy) : x(xx), y(yy) {}

	vec2 operator*(const double& a) const {
		return vec2(x*a, y*a);
	}
	vec2 operator*(const vec2& a) const {
		return vec2(a.x*x, a.y*y);
	}
	vec2 operator/(const vec2& a) const {
		return vec2(a.x/x, a.y/y);
	}
	vec2 operator+(const vec2& a) const {
		return vec2(a.x+x, a.y+y);
	}
	vec2 operator-(const vec2& a) const {
		return vec2(a.x-x, a.y-y);
	}
};

double dist(vec2 a, vec2 b) {
	return sqrt(((b.x - a.x) * (b.x - a.x)) + ((b.y - a.y) * (b.y - a.y)));
}

struct color {
	int r, g, b, a;
	Uint32 c;

	color() {}
	color(Uint32 col) : c(col) {}
	color(int red, int green, int blue, int alpha) : r(red), g(green), b(blue), a(alpha) {
		c = (alpha << 24) | (red << 16) | (green << 8) | blue;
	}

	Uint32 darken() {
		double m = 0.60;
		return (int(a * m) << 24) | (int(r * m) << 16) | (int(g * m) << 8) | int(b * m);
	}
};

struct timeStep {
	int now, last;
	double time;

	timeStep(int n) : now(n), last(0), time(0) {}

	void update(int n, int freq) {
		last = now;
		now = n;
		time = (double)((now - last) / (double)freq);
	}
};

struct camera {
	vec2 position;
	vec2 direction;
	vec2 plane;
	double moveSpeed, turnSpeed;
	bool left, right, forward, backward;

	camera(vec2 pos, vec2 dir, vec2 pln, double ms, double ts) : position(pos), direction(dir), plane(pln), moveSpeed(ms), turnSpeed(ts) {}
};

struct texture {
	std::string filename;
	int width, height;
	int id;
	bool transparent;
	std::vector<color> data;

	texture() {}
	texture(std::string fn, int id, bool ts) : filename(fn), id(id), transparent(ts) {
	    SDL_Surface *surface = SDL_LoadBMP(filename.c_str());

	    width = surface->w;
	    height = surface->h;

	    SDL_Surface *rgbaSurface = SDL_ConvertSurfaceFormat(surface, SDL_PIXELFORMAT_RGBA32, 0);
	    SDL_FreeSurface(surface);
	    surface = rgbaSurface;

	    Uint32 *pixels = static_cast<Uint32*>(surface->pixels);
	    data.reserve(width * height);

	    for (int y = 0; y < height; ++y) {
	        for (int x = 0; x < width; ++x) {
	            Uint32 pixelColor = pixels[y * width + x];
	            Uint8 r, g, b, a;
	            SDL_GetRGBA(pixelColor, surface->format, &r, &g, &b, &a);
	            data.emplace_back(r, g, b, a);
	        }
	    }

	    SDL_FreeSurface(surface);
}
};

struct ray {
	vec2 position; 		// Initial start position of the ray in floating point value
	vec2 tilePosition; 	// Initial start position of the ray in tile world space

	vec2 direction;		// Direction that the ray shoots
	
	vec2 initialStepLength; // The initial length a ray has to travel to reach the edge of its current tile
	vec2 stepDirection; 	// The x or y direction the ray will step in to the next tile
	vec2 sideLength; 		// The side lengths of the ray traveling through a tile

	vec2 intersectionPoint; // The point at where the ray intersects an object
	vec2 intersectedTile;

	int currentTravelLength, maxTravelLength; // making a maximum ray travel distance
	double length; // length of ray from starting position to its intersection point
	double distance; // distance of the ray from its starting position to the intersected tile
	bool instersectedSide; // true = right side, false = left side
	int intersectedNumber; // the number that the ray intersected in the array

	int lineStart, lineEnd, lineHeight;

	std::vector<ray> intersections;
	bool multiIntersections;

	ray() {}
	ray(camera& cam, double camSpcX) : currentTravelLength(0), maxTravelLength(100) {
		position = cam.position;
		tilePosition = vec2(int(cam.position.x), int(cam.position.y));
		direction = cam.direction + cam.plane * camSpcX;
		initialStepLength = vec2((direction.x == 0) ? 1e30 : std::abs(1 / direction.x), (direction.y == 0) ? 1e30 : std::abs(1 / direction.y));

		if (direction.x < 0) {
			stepDirection.x = -1;
			sideLength.x = (position.x - tilePosition.x) * initialStepLength.x;
		}
		else {
			stepDirection.x = 1;
			sideLength.x = (tilePosition.x + 1.0 - position.x) * initialStepLength.x;
		}

		if (direction.y < 0) {
			stepDirection.y = -1;
			sideLength.y = (position.y - tilePosition.y) * initialStepLength.y;
		}
		else {
			stepDirection.y = 1;
			sideLength.y = (tilePosition.y + 1.0 - position.y) * initialStepLength.y;
		}
	}

	void findIntersection() {
		while (true) {
			// checking to see if ray has exceeded its max length
			if (currentTravelLength > maxTravelLength) {
				break;
			}

			// jumping to next tile in either x or y direction
			if (sideLength.x < sideLength.y) {
				tilePosition.x += stepDirection.x;
				sideLength.x += initialStepLength.x;
				instersectedSide = true;
			}
			else {
				tilePosition.y += stepDirection.y;
				sideLength.y += initialStepLength.y;
				instersectedSide = false;
			}

			if (m[int(tilePosition.x)][int(tilePosition.y)] > 0) {
				intersectedNumber = m[int(tilePosition.x)][int(tilePosition.y)];
				intersectedTile = vec2(int(tilePosition.x), int(tilePosition.y));

				if (instersectedSide) {
					length = sideLength.x - initialStepLength.x;
				}
				else {
					length = sideLength.y - initialStepLength.y;
				}

				intersectionPoint.x = position.x + length * direction.x;
				intersectionPoint.y = position.y + length * direction.y;

				distance = dist(position, intersectionPoint);

				if (intersectedNumber == 5 || intersectedNumber == 6) {
					intersections.push_back(*this);
					multiIntersections = true;
					continue;
				}

				break;
			}
			currentTravelLength += 1;
		}

		lineHeight = (HEIGHT / length);

		// calculate lowest and highest pixel to fill in current stripe
		lineStart = -lineHeight / 2 + HEIGHT / 2;
		lineEnd = lineHeight / 2 + HEIGHT / 2;

		if (lineStart < 0) {
			lineStart = 0;
		}

		if (lineEnd >= HEIGHT) {
			lineEnd = HEIGHT - 1;
		}
	}
};

struct world {
	int map[24][24];

	std::vector<texture> textures; // needs to be in order from 0

	texture currentTexture;

	vec2 hitDistance;
	vec2 texCoord;

	double step;
	double texPos;

	world(int map[][24], std::vector<texture> textures) : map{{0}}, textures(textures) {
		for (int i = 0; i < 24; i++) {
			for (int j = 0; j < 24; j++) {
				this->map[i][j] = map[i][j];
			}
		}
	}

	void findCollision(ray r) {
		vec2 hitDistance = r.intersectionPoint - r.intersectedTile;
		for (int i = 0; i < textures.size(); ++i) {
			if (map[int(r.tilePosition.x)][int(r.tilePosition.y)] == textures[i].id) {
				if (r.instersectedSide) {
					currentTexture = textures[i];
					texCoord.x = std::abs(hitDistance.y) * currentTexture.width;
				}
				else {
					currentTexture = textures[i];
					for (int j = 0; j < currentTexture.data.size(); ++j) {
						currentTexture.data[j].darken();
					}
					texCoord.x = std::abs(hitDistance.x) * currentTexture.width;
				}
			}
		}

		step = 1.0 * currentTexture.width / r.lineHeight;
		// Starting texture coordinate
		texPos = (r.lineStart - HEIGHT / 2 + r.lineHeight / 2) * step;
	}

	void updateY() {
		texCoord.y = (int)texPos & (currentTexture.width - 1);
		texPos += step;
	}
};

int main() {
	SDL_Init(SDL_INIT_VIDEO);

	SDL_Window *win = SDL_CreateWindow("Renderer", 100, 100, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
	SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

	timeStep ts = timeStep(SDL_GetPerformanceCounter());
	float time = 0.0f;

	camera player = camera(vec2(2.0, 2.0), vec2(-1, 0), vec2(0, 0.66), 5, 3);

	std::vector<texture> textures;

	texture dirt = 					texture(std::string("res/dirt.bmp"), 			0, false);
	texture plank = 				texture(std::string("res/oak_planks.bmp"),	1, false);
	texture cobblestone = 	texture(std::string("res/cobblestone.bmp"),	2, false);
	texture log = 					texture(std::string("res/oak_log.bmp"), 		3, false);
	texture brick = 				texture(std::string("res/brick.bmp"), 		4, false);
	texture normal_glass =	texture(std::string("res/glass.bmp"), 		5, true);

	textures.push_back(plank);
	textures.push_back(cobblestone);
	textures.push_back(dirt);
	textures.push_back(log);
	textures.push_back(brick);

	textures.push_back(normal_glass);

	world map = world(m, textures);

	SDL_Event e;
	bool quit = false;
	while (!quit){
		ts.update(SDL_GetPerformanceCounter(), SDL_GetPerformanceFrequency());

		time += ts.time;
		if (time > 1.0f) {
			SDL_SetWindowTitle(win, std::to_string(1.0 / ts.time).c_str());
			time = 0.0f;
		}
		// Create a surface from the bitmap array
		SDL_Surface *surface = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, 0, 0, 0, 0);

		// make top and bottom colors
		color bottom = color(113, 113, 113, 255);
		color top = color(56, 56, 56, 255);
		Uint32 *pixels = (Uint32 *)surface->pixels;
		for (int y = 0; y < HEIGHT; y++) {
			for (int x = 0; x < WIDTH; x++) {
				if (y * WIDTH + x > (WIDTH * HEIGHT) / 2) {
					pixels[y * WIDTH + x] = top.c;
				}
				else {
					pixels[y * WIDTH + x] = bottom.c;
				}
			}
		}


		for (int x = 0; x < WIDTH; ++x) {
			double cameraSpaceX = ((2 * x) / double(WIDTH)) - 1; // camera x coordinate ranging from -1(left) to 1(right)

			ray r = ray(player, cameraSpaceX);
			r.findIntersection();

			map.findCollision(r);

			for (int y = 0; y < HEIGHT; ++y) {
				if (y > r.lineStart && y < r.lineEnd) {
					map.updateY();

					if (r.instersectedSide) {
						pixels[y * WIDTH + x] = map.currentTexture.data[map.texCoord.y * map.currentTexture.width + map.texCoord.x].darken();
					}
					else {
						pixels[y * WIDTH + x] = map.currentTexture.data[map.texCoord.y * map.currentTexture.width + map.texCoord.x].c;
					}
				}
			}

			if (r.intersections.size() > 0) {
				for (int i = 0; i < r.intersections.size(); ++i) {
					ray i_r = r.intersections[i];

					world i_map = world(m, textures);
					i_map.findCollision(i_r);

					int newWallHeight = (HEIGHT / i_r.length);

					int newLineStart = -newWallHeight / 2 + HEIGHT / 2;
					int newLineEnd = newWallHeight / 2 + HEIGHT / 2;

					if (newLineStart < 0) {
						newLineStart = 0;
					}
					if (newLineEnd >= HEIGHT) {
						newLineEnd = HEIGHT - 1;
					}

					vec2 texCoord = i_r.intersectionPoint - i_r.intersectedTile; // if instersectedSide then use y, else use x
					int xTex;

					if (i_r.instersectedSide) {
						xTex = double(std::abs(texCoord.y)) * double(i_map.currentTexture.width);
					}
					else {
						xTex = double(std::abs(texCoord.x)) * double(i_map.currentTexture.width);
					}

					double step = 1.0 * i_map.currentTexture.width / newWallHeight;
					// Starting texture coordinate
					double texPos = (newLineStart - HEIGHT / 2 + newWallHeight / 2) * step;


					for (int y = 0; y < HEIGHT; ++y) {
						if (y > newLineStart && y < newLineEnd) {
							int proportion = ((double(y) - double(newLineStart)) / (double(newLineEnd) - double(newLineStart))) * i_map.currentTexture.width;

							int texY = (int)texPos & (16 - 1);
							texPos += step;

							i_map.updateY();

							if (i_map.currentTexture.data[texY * i_map.currentTexture.width + xTex].a != 0) {
								pixels[y * WIDTH + x] = i_map.currentTexture.data[texY * i_map.currentTexture.width + xTex].c;
							}	
						}
					}
				}
			}
		}

		SDL_Texture* texture = SDL_CreateTextureFromSurface(ren, surface);

		if (player.forward) {
			if(m[int(player.position.x + player.direction.x * (player.moveSpeed * ts.time))][int(player.position.y)] == false) player.position.x += player.direction.x * (player.moveSpeed * ts.time);
			if(m[int(player.position.x)][int(player.position.y + player.direction.y * (player.moveSpeed * ts.time))] == false) player.position.y += player.direction.y * (player.moveSpeed * ts.time);
		}

		if (player.backward) {
			if(m[int(player.position.x - player.direction.x * (player.moveSpeed * ts.time))][int(player.position.y)] == false) player.position.x -= player.direction.x * (player.moveSpeed * ts.time);
			if(m[int(player.position.x)][int(player.position.y - player.direction.y * (player.moveSpeed * ts.time))] == false) player.position.y -= player.direction.y * (player.moveSpeed * ts.time);
		}

		if (player.left) {
			vec2 oldDirection = player.direction;
			vec2 oldPlane = player.plane;

			// rotate left
			player.direction.x = player.direction.x * cos(player.turnSpeed * ts.time) - player.direction.y * sin(player.turnSpeed * ts.time);
			player.direction.y = oldDirection.x * sin(player.turnSpeed * ts.time) + player.direction.y * cos(player.turnSpeed * ts.time);
			player.plane.x = player.plane.x * cos(player.turnSpeed * ts.time) - player.plane.y * sin(player.turnSpeed * ts.time);
			player.plane.y = oldPlane.x * sin(player.turnSpeed * ts.time) + player.plane.y * cos(player.turnSpeed * ts.time);
		}

		if (player.right) {
			vec2 oldDirection = player.direction;
			vec2 oldPlane = player.plane;

			// rotate the right
			player.direction.x = player.direction.x * cos(-player.turnSpeed * ts.time) - player.direction.y * sin(-player.turnSpeed * ts.time);
			player.direction.y = oldDirection.x * sin(-player.turnSpeed * ts.time) + player.direction.y * cos(-player.turnSpeed * ts.time);
			player.plane.x = player.plane.x * cos(-player.turnSpeed * ts.time) - player.plane.y * sin(-player.turnSpeed * ts.time);
			player.plane.y = oldPlane.x * sin(-player.turnSpeed * ts.time) + player.plane.y * cos(-player.turnSpeed * ts.time);
		}

		SDL_RenderClear(ren);
		//Draw the texture
		SDL_RenderCopy(ren, texture, NULL, NULL);
		//Update the screen
		SDL_RenderPresent(ren);

		SDL_UpdateWindowSurface(win); 

		while (SDL_PollEvent(&e)){
			if (e.type == SDL_QUIT){
				quit = true;
			}
			if (e.type == SDL_KEYDOWN){
				if (e.key.keysym.sym == SDLK_ESCAPE) {
					quit = true;
				}
				if (e.key.keysym.sym == SDLK_w) {
					player.forward = true;
				} 
				if (e.key.keysym.sym == SDLK_s) {
					player.backward = true;
				} 
				if (e.key.keysym.sym == SDLK_a) {
					player.left = true;
				} 
				if (e.key.keysym.sym == SDLK_d) {
					player.right = true;
				} 

			}
			if (e.type == SDL_KEYUP){
				if (e.key.keysym.sym == SDLK_w) {
					player.forward = false;
				} 
				if (e.key.keysym.sym == SDLK_s) {
					player.backward = false;
				} 
				if (e.key.keysym.sym == SDLK_a) {
					player.left = false;
				} 
				if (e.key.keysym.sym == SDLK_d) {
					player.right = false;
				} 
			}
		}
		SDL_FreeSurface(surface);
		SDL_DestroyTexture(texture);
	}
	SDL_DestroyRenderer(ren);
	SDL_DestroyWindow(win);
	SDL_Quit();

	return 0;
}

// uint32_t pixel;
// uint8_t r = 0;
// uint8_t g = 0; 
// uint8_t b = 0; 
// uint8_t a = 255;
// pixel = (a << 24) | (r << 16) | (g << 8) | b;