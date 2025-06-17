# ~~~~
'''
Python program to catalog legendrian knot mosaics from a file based on Thurston Bennequin, rotation, and smooth knot type invariants.
Mosaics are represented by a base 10 number, which is read from the file and converted into a list.
The program then traverses the knot from a starting tile, calculating tb/rotation along the way.
If every tile has been traversed by the time we return to our starting point, we know that the mosaic represents a knot.

We also build the extended gauss code of the knot during traversal, which is used to calculate the knot's HOMFLY as detailed in the sagemath link below.
The invariants of this new knot are then compared to a list of knots found previously, and it is added to this list if its combination of invatiants have not occurred before.

The total number of mosaics representing knots is printed when the program ends, and the knot catalog is written to a file.

Note: Every input file should start with an empty mosaic of the same size as the other mosaics, e.g:
000000000
000021034
000210340
...
Files produced by the associated rust program will always start this way.

Note: This program should be run in sagemath, which can be done by calling the program with sage from a terminal as one would with python (sage file_cat.py)

sagemath's Links library is used to calculate homfly polynomials -- see (https://doc.sagemath.org/html/en/reference/knots/sage/knots/link.html) for more details
See the associated rust file for information on how these mosaics are generated.

This material is based upon work supported by the National Science Foundation under Grant No. MPS-2150299
'''

from sage.all import *

#Basic driver function
def main():
    print("file to read from:")
    in_file = input()
    print('output file?:')
    out_file = input()

    cylindrical_mosaic.file_catalog(in_file, out_file, False)

class cylindrical_mosaic:
    #Valid connections for each tile
    #The first digit in the tuple represents the incoming face, the second the outgoing face.
    #Faces are assigned as below:
    #    1
    # 2 ▇▇ 0
    #    3
    valid_connections = (
        (()),
        ((2,3),(3,2)),
        ((0,3),(3,0)),
        ((0,1),(1,0)),
        ((2,1),(1,2)),
        ((2,0),(0,2)),
        ((1,3),(3,1)),
        ((2,3),(3,2),(1,0),(0,1)),
        ((0,3),(3,0),(2,1),(1,2)),
        ((0,2),(1,3),(2,0),(3,1)),
        ((0,2),(1,3),(2,0),(3,1)),
        ((0,2),(1,3),(2,0),(3,1))
    )

    @classmethod
    def file_catalog(cls, input_name, output_name, images):
        knot_catalog = dict()

        input_file = open(input_name, 'r')
        output_file = open(output_name, 'w')

        test_string = input_file.readline().strip()
        size = int(len(test_string)**(0.5))
        mosaic = [10]*(size**2)
        satisfied = [False]*(size ** 2)
        crossing_strands = [[0]*4 for _ in range((size ** 2))]
        made_connections = [[] for _ in range((size ** 2))]
        crossing_indices = []
        pd_codes = []
        curr_tile = 0
        starting_tile = None
        face = 0
        strand_number = 1
        i = 0
        num = 0
        knot_count = 0

        knot = None
        for mosaic_string in input_file:
            made_connections = [[] for _ in range(size ** 2)]
            crossing_indices = []
            pd_codes = []
            curr_tile = int(0)
            starting_tile = None
            strand_number = 1

            k = num = 0
            for char in mosaic_string.strip():
                num = int(char, base = 16)
                mosaic[k] = num
                satisfied[k] = num == 0
                if starting_tile == None and num != 0:
                    starting_tile = k
                k += 1

            curr_tile = starting_tile
            face = cls.valid_connections[mosaic[curr_tile]][0][0]
            not_looped = True
            while not_looped:
                for conn in cls.valid_connections[mosaic[curr_tile]]:
                    if conn[0] == face:
                        if conn in made_connections[curr_tile]:
                            not_looped = False
                            break
                        made_connections[curr_tile].append(conn)
                        if ((len(made_connections[curr_tile]) == 1) and mosaic[curr_tile] < 7) or (len(made_connections[curr_tile]) == 2):
                            satisfied[curr_tile] = True

                        #Crossing logic
                        if mosaic[curr_tile] > 8:
                            if satisfied[curr_tile]:
                                crossing_indices.append(curr_tile)
                            crossing_strands[curr_tile][face] = strand_number
                            strand_number += 1
                            crossing_strands[curr_tile][(face + 2) % 4] = strand_number
                        else:
                            face = (conn[1] + 2) % 4

                        #Go to next tile
                        if face == 0: #left
                            curr_tile = size*((curr_tile)//size) + (curr_tile + size - 1)%size
                        elif face == 1: #down
                            curr_tile += size
                        elif face == 2: #right
                            curr_tile = size*((curr_tile)//size) + (curr_tile + 1)%size
                        elif face == 3: #up
                            curr_tile -= size
                        break

            if all(satisfied):
                knot_count += 1
                if len(crossing_indices) < 3 :
                    homf = 1
                else:
                    for index in crossing_indices:
                        if mosaic[index] == 9:
                            if (0,2) in made_connections[index]:
                                pd_codes.append(crossing_strands[index])
                            else:
                                pd_codes.append(crossing_strands[index][2:] + crossing_strands[index][:2]) #Rotated by 2
                        else:
                            if (1,3) in made_connections[index]:
                                pd_codes.append(crossing_strands[index][1:] + crossing_strands[index][:1])
                            else:
                                pd_codes.append(crossing_strands[index][3:] + crossing_strands[index][:3])
                    for l in range(len(pd_codes[-1])):
                        if pd_codes[-1][l] == strand_number:
                            pd_codes[-1][l] = 1
                            break
                    homf = Link(pd_codes).homfly_polynomial()
                if not homf in knot_catalog:
                    knot_catalog[homf] = mosaic_string.strip()
                    output_file.write(f"\t{homf}: {knot_catalog[homf]}\n")
                    output_file.flush()
                    if images:
                        to_png(mosaic,f"images/{output_name}/{mosaic_string.strip()}.png")
        print(knot_count)
        output_file.close()

def to_png(matrix,output_filename):
    tile_size = 64
    border_size = 4
    border_color = (196, 196, 196, 255)
    size = int(len(matrix)**0.5)

    tile_images = {}
    for num in range(11):
        file_name = f"tiles/{num}.png"
        try:
            tile_images[num] = Image.open(file_name).convert("RGBA")
        except FileNotFoundError:
            print(f"Failed to load image {file_name}")

    mosaic_width = size * tile_size + 2 * border_size
    mosaic = Image.new("RGBA", (mosaic_width, mosaic_width), border_color)
    draw = ImageDraw.Draw(mosaic)

    for i, tile in enumerate(matrix):
            if tile in tile_images:
                img_tile = tile_images[tile]
                for y in range(tile_size):
                    for x in range(tile_size):
                        pixel = img_tile.getpixel((x, y))
                        mosaic.putpixel(( (i % size) * tile_size + x + border_size, (i // size) * tile_size + y + border_size), pixel)
    mosaic.save(output_filename)
main()
