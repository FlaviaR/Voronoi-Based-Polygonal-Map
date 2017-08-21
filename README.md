# Voronoi-Based-Polygonal-Map

This project uses Voronoi diagrams for polygonal map generation. 


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This project uses CGAL and OpenGL.


## How to Run 

* Open a terminal and go to the downloaded directory
* This project includes a makefile and can be run on different operating systems: [linux, osx (mac), win32 (32-bit windows)]
* The default operating system is linux, so if that is the OS that you are running, then simply type 'make'
* If not, then type 'make OS=[desired OS]'. For example, if using a mac, type 'make OS=osx'
* Run the generated executable file

## Keyboard and Mouse Interactivity
* 'h' - display help
* '1' - Initial Voronoi diagram
* '2' - Voronoi diagram with applied Lloyd Relaxation
* '3' - Radially generated map
* 'r' - randomize - generate a new map
* 'MOUSE LEFT BUTTON' - display Delaunay triangulation
* 'MOUSE RIGHT BUTTON' - display Voronoi diagram


## TO DO:
* Finish developing graph structure
* Display lakes
* Display rivers
* Add elevations
* Add biomes
* Add town names
* Add other map generation algs 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Inspiration - "Polygonal Map Generation for Games" from Red Blob Games.
* See - http://www-cs-students.stanford.edu/~amitp/game-programming/polygon-map-generation/
