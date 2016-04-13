#Shredder

Shredder is a command line tool that takes a 3D mesh as input and produces a series of SVG files with laser-cut-ready contours. You can either indicate the number of slices in which you want to cut your 3D mesh -- shredder will the compute the slice thickness, or indicate the slice thickness and let shredder compute the number of slices. You can also indicate the maximum width and height of each of the SVG files. Based on those sizes, shredder will compute the number of slices that can be fit in each SVG file.

The command lines look like this:

````
shredder
	input_mesh.txt      // input text format mesh
    output_root         // name root to be used for the svg files
    -n number_of_slices // number of slices
    -w maximum_width_in_mm
    -h maximum_height_in_mm
    -s scale
````

````
shredder
	input_mesh.txt                // input text format mesh
    output_root                   // name root to be used for the svg files
    -t thickness_of_slices_in_mm  // thickness of slices
    -w maximum_width_in_mm
    -h maximum_height_in_mm
    -s scale
````

#Example
As an example, we will slice the brain of a cat, obtained from the Brain Catalogue (http://braincatalogue.org/Cat). The surface is in .ply format. The first step is to convert it to the .txt format that shredder can read. For this we use meshgeometry (http://github.com/r03ert0/meshgeometry):


````
meshgeometry_mac -i example/cat.ply -o example/cat.txt
````

Now, we use shredder to slice this brain in 50 slices (the slice thickness will be automatically computed). The slices will be layed out in pages of width x height = 300 x 200:

````
shredder example/cat.txt cat -n 50 -w 300 -h 200 -s 1
````

This command produces 49 individual slices, layed out in 2 files. Each slice is numbered. You can use toothpicks to stack the slices back together using the grid of small square holes.